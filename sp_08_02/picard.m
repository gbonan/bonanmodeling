function [soil] = picard (soil, params, dt, dpsi_tolerance, water_balance_error)

% -------------------------------------------------------------
% Use modified Picard iteration to solve the Richards equation
% -------------------------------------------------------------

% Input
% dt                   ! Time step (s)
% dpsi_tolerance       ! Convergence criterion for dpsi
% water_balance_error  ! Water balance error tolerance
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

% --- Initialization

for i = 1:soil.nsoi

   % Save current soil moisture for time n
   theta0(i) = soil.theta(i);

   % Initialize delta_psi to a large value
   dpsi(i) = 1.0e36;

end

% --- Iteration loop

m = 0;
while (max(abs(dpsi)) > dpsi_tolerance) 

   % Increment iteration counter

   m = m + 1;

   % Stop if too many iterations

   if (m > 50)
      error ('Too many iterations')
   end

   % Hydraulic properties for current psi
   % theta - volumetric soil moisture
   % K -  hydraulic conductivity
   % cap - specific moisture capacity

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
   b(i) = soil.cap(i) * soil.dz(i) / dt - a(i) - c(i);
   d(i) = K_onehalf / dz_onehalf * (soil.psi0 - soil.psi(i)) ... 
        - K_plus_onehalf(i) / soil.dz_plus_onehalf(i) * (soil.psi(i) - soil.psi(i+1)) ...
        + K_onehalf - K_plus_onehalf(i) - (soil.theta(i) - theta0(i)) * soil.dz(i) / dt;

   for i = 2:soil.nsoi-1
      a(i) = -K_plus_onehalf(i-1) / soil.dz_plus_onehalf(i-1);
      c(i) = -K_plus_onehalf(i) / soil.dz_plus_onehalf(i);
      b(i) = soil.cap(i) * soil.dz(i) / dt - a(i) - c(i);
      d(i) = K_plus_onehalf(i-1) / soil.dz_plus_onehalf(i-1) * (soil.psi(i-1) - soil.psi(i)) ...
           - K_plus_onehalf(i) / soil.dz_plus_onehalf(i) * (soil.psi(i) - soil.psi(i+1)) ...
           + K_plus_onehalf(i-1) - K_plus_onehalf(i) - (soil.theta(i) - theta0(i)) * soil.dz(i) / dt;
   end

   i = soil.nsoi;
   a(i) = -K_plus_onehalf(i-1) / soil.dz_plus_onehalf(i-1);
   c(i) = 0;
   b(i) = soil.cap(i) * soil.dz(i) / dt - a(i) - c(i);
   d(i) = K_plus_onehalf(i-1) / soil.dz_plus_onehalf(i-1) * (soil.psi(i-1) - soil.psi(i)) ...
        + K_plus_onehalf(i-1) - soil.K(i) - (soil.theta(i) - theta0(i)) * soil.dz(i) / dt;

   % Solve for the change in psi

   [dpsi] = tridiagonal_solver (a, b, c, d, soil.nsoi);

   % Update psi

   for i = 1:soil.nsoi
      soil.psi(i) = soil.psi(i) + dpsi(i);
   end

end

% --- Check water balance

soil.Q0 = -K_onehalf / dz_onehalf * (soil.psi0 - soil.psi(1)) - K_onehalf;
soil.QN = -soil.K(soil.nsoi);

soil.dtheta = 0;
for i = 1:soil.nsoi
   soil.dtheta = soil.dtheta + (soil.theta(i) - theta0(i)) * soil.dz(i);
end

soil.err = soil.dtheta - (soil.QN - soil.Q0) * dt;
if (abs(soil.err) > water_balance_error)
   error ('Water conservation error')
end
