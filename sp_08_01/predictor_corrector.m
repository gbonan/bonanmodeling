function [soil] = predictor_corrector (soil, params, dt)

% -------------------------------------------------------------
% Use predictor-corrector method to solve the Richards equation
% -------------------------------------------------------------

% Input
% dt                   ! Time step (s)
% soil.nsoi            ! Number of soil layers
% soil.functions       ! van Genuchten or Campbell relationships
% soil.dz_plus_onehalf ! Thickness between between z(i) and z(i+1) (cm)
% soil.dz              ! Soil layer thickness (cm)
% soil.psi0            ! Soil surface matric potential boundary condition (cm)
%
% Input/output
% soil.theta           ! Volumetric soil moisture
% soil.psi             ! Matric potential (cm)
%
% Output
% soil.K               ! Hydraulic conductivity (cm H2O/s)
% soil.cap             ! Specific moisture capacity (/cm)
% soil.Q0              ! Infiltration flux (cm H2O/s)
% soil.QN              ! Drainage flux (cm H2O/s)
% soil.dtheta          ! Change in soil moisture (cm H2O)
% soil.err             ! Water balance error (cm H2O)

% --- Save current soil moisture and matric potential for time n

for i = 1:soil.nsoi
   theta_n(i) = soil.theta(i);
   psi_n(i) = soil.psi(i);
end

% --- Predictor step using implict solution for time n+1/2

% Hydraulic properties for current psi:
% theta - volumetric soil moisture
% K     - hydraulic conductivity
% cap   - specific moisture capacity

for i = 1:soil.nsoi
   switch soil.functions
      case 'van_Genuchten'
      [soil.theta(i), soil.K(i), soil.cap(i)] = van_Genuchten (params, soil.psi(i));
      case 'Campbell'
      [soil.theta(i), soil.K(i), soil.cap(i)] = Campbell (params, soil.psi(i));
    end
end

% Hydraulic conductivity at i+1/2 interface between layers i and i+1 is the arithmetic mean

for i = 1:soil.nsoi-1
   K_plus_onehalf(i) = 0.5 * (soil.K(i) + soil.K(i+1));
end

% Hydraulic conductivity at i=1/2 between surface (i=0) and first layer i=1

K_onehalf = soil.K(1);

% dz at i=1/2 between surface (i=0) and first layer i=1

dz_onehalf = 0.5 * soil.dz(1);

% Terms for tridiagonal matrix

i = 1;
a(i) = 0;
c(i) = -K_plus_onehalf(i) / soil.dz_plus_onehalf(i);
b(i) = soil.cap(i) * soil.dz(i) / (0.5 * dt) + K_onehalf / dz_onehalf - c(i);
d(i) = soil.cap(i) * soil.dz(i) / (0.5 * dt) * soil.psi(i) + K_onehalf / dz_onehalf * soil.psi0 + K_onehalf - K_plus_onehalf(i);

for i = 2:soil.nsoi-1
   a(i) = -K_plus_onehalf(i-1) / soil.dz_plus_onehalf(i-1);
   c(i) = -K_plus_onehalf(i) / soil.dz_plus_onehalf(i);
   b(i) = soil.cap(i) * soil.dz(i) / (0.5 * dt) - a(i) - c(i);
   d(i) = soil.cap(i) * soil.dz(i) / (0.5 * dt) * soil.psi(i) + K_plus_onehalf(i-1) - K_plus_onehalf(i);
end

i = soil.nsoi;
a(i) = -K_plus_onehalf(i-1) / soil.dz_plus_onehalf(i-1);
c(i) = 0;
b(i) = soil.cap(i) * soil.dz(i) / (0.5 * dt) - a(i) - c(i);
d(i) = soil.cap(i) * soil.dz(i) / (0.5 * dt) * soil.psi(i) + K_plus_onehalf(i-1) - soil.K(i);

% Solve for psi at n+1/2

[psi_pred] = tridiagonal_solver (a, b, c, d, soil.nsoi);

% --- Corrector step using Crank-Nicolson solution for time n+1

% Hydraulic properties for psi_pred

for i = 1:soil.nsoi
   switch soil.functions
      case 'van_Genuchten'
      [soil.theta(i), soil.K(i), soil.cap(i)] = van_Genuchten (params, psi_pred(i));
      case 'Campbell'
      [soil.theta(i), soil.K(i), soil.cap(i)] = Campbell (params, psi_pred(i));
    end
end

% Hydraulic conductivity at i+1/2 interface between layers i and i+1

for i = 1:soil.nsoi-1
   K_plus_onehalf(i) = 0.5 * (soil.K(i) + soil.K(i+1));
end

% Hydraulic conductivity at i=1/2 between surface (i=0) and first layer i=1

K_onehalf = soil.K(1);

% dz at i=1/2 between surface (i=0) and first layer i=1

dz_onehalf = 0.5 * soil.dz(1);

% Terms for tridiagonal matrix

i = 1;
a(i) = 0;
c(i) = -K_plus_onehalf(i) / (2 * soil.dz_plus_onehalf(i));
b(i) = soil.cap(i) * soil.dz(i) / dt  + K_onehalf / (2 * dz_onehalf) - c(i);
d(i) = soil.cap(i) * soil.dz(i) / dt * soil.psi(i) + K_onehalf / (2 * dz_onehalf) * soil.psi0 ...
     + K_onehalf / (2 * dz_onehalf) * (soil.psi0 - soil.psi(i)) ...
     + c(i) * (soil.psi(i) - soil.psi(i+1)) + K_onehalf - K_plus_onehalf(i);

for i = 2:soil.nsoi-1
   a(i) = -K_plus_onehalf(i-1) / (2 * soil.dz_plus_onehalf(i-1));
   c(i) = -K_plus_onehalf(i) / (2 * soil.dz_plus_onehalf(i));
   b(i) = soil.cap(i) * soil.dz(i) / dt - a(i) - c(i);
   d(i) = soil.cap(i) * soil.dz(i) / dt * soil.psi(i) - a(i) * (soil.psi(i-1) - soil.psi(i)) ...
        + c(i) * (soil.psi(i) - soil.psi(i+1)) + K_plus_onehalf(i-1) - K_plus_onehalf(i);
end

i = soil.nsoi;
a(i) = -K_plus_onehalf(i-1) / (2 * soil.dz_plus_onehalf(i-1));
c(i) = 0;
b(i) = soil.cap(i) * soil.dz(i) / dt - a(i) - c(i);
d(i) = soil.cap(i) * soil.dz(i) / dt * soil.psi(i) - a(i) * (soil.psi(i-1) - soil.psi(i)) + K_plus_onehalf(i-1) - soil.K(i);

% Solve for psi at n+1

[soil.psi] = tridiagonal_solver (a, b, c, d, soil.nsoi);

% --- Check water balance

soil.Q0 = -K_onehalf / (2 * dz_onehalf) * ((soil.psi0 - psi_n(1)) + (soil.psi0 - soil.psi(1))) - K_onehalf;
soil.QN = -soil.K(soil.nsoi);

soil.dtheta = 0;
for i = 1:soil.nsoi
   soil.dtheta = soil.dtheta + (soil.theta(i) - theta_n(i)) * soil.dz(i);
end

soil.err = soil.dtheta - (soil.QN - soil.Q0) * dt;
