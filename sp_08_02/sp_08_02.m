% Supplemental program 8.2

% ---------------------------------------------------------------------
% Use modified Picard iteration to solve the Richards equation for
% infiltration with surface soil moisture as the boundary condition.
% ---------------------------------------------------------------------

% --- Define soil layers

% Number of soil layers

soil.nsoi = 150;

% Soil layer thickness (cm)

for i = 1:soil.nsoi
   soil.dz(i) = 1.0;
end

% Soil depth (cm) at i+1/2 interface between layers i and i+1 (negative distance from surface)

soil.z_plus_onehalf(1) = -soil.dz(1);
for i = 2:soil.nsoi
   soil.z_plus_onehalf(i) = soil.z_plus_onehalf(i-1) - soil.dz(i);
end

% Soil depth (cm) at center of layer i (negative distance from surface)

soil.z(1) = 0.5 * soil.z_plus_onehalf(1);
for i = 2:soil.nsoi
   soil.z(i) = 0.5 * (soil.z_plus_onehalf(i-1) + soil.z_plus_onehalf(i));
end

% Thickness between between z(i) and z(i+1)

for i = 1:soil.nsoi-1
   soil.dz_plus_onehalf(i) = soil.z(i) - soil.z(i+1);
end
soil.dz_plus_onehalf(soil.nsoi) = 0.5 * soil.dz(soil.nsoi);

% --- Soil parameters

 soil.functions = 'van_Genuchten';  % Use van Genuchten relationships
%soil.functions = 'Campbell';       % Use Campbell relationships

switch soil.functions

   case 'Campbell'
   % example from Hornberger & Wiberg (2005, Fig. 8.3)
   ityp = 0;              % Soil texture flag
   theta_sat = 0.25;      % Volumetric water content at saturation
   psi_sat = -25.0;       % Matric potential at saturation (cm)
   bc = 0.2;              % Exponent
   Ksat = 3.4e-03;        % Hydraulic conductivity at saturation (cm/s)
   params = [theta_sat psi_sat bc Ksat ityp];

   case 'van_Genuchten'

   % Haverkamp et al. (1977): sand !!! DOES NOT WORK !!!
%  ityp = 1;              % Soil texture flag
%  theta_res = 0.075;     % Residual water content
%  theta_sat = 0.287;     % Volumetric water content at saturation
%  vg_alpha = 0.027;      % Inverse of the air entry potential (/cm)
%  vg_n = 3.96;           % Pore-size distribution index
%  vg_m = 1;              % Exponent
%  Ksat = 34 / 3600;      % Hydraulic conductivity at saturation (cm/s)

   % Haverkamp et al. (1977): Yolo light clay
   ityp = 2;              % Soil texture flag
   theta_res = 0.124;     % Residual water content
   theta_sat = 0.495;     % Volumetric water content at saturation
   vg_alpha = 0.026;      % Inverse of the air entry potential (/cm)
   vg_n = 1.43;           % Pore-size distribution index
   vg_m = 1 - 1 / vg_n;   % Exponent
   Ksat = 0.0443 / 3600;  % Hydraulic conductivity at saturation (cm/s)

   params = [theta_res theta_sat vg_alpha vg_n vg_m Ksat ityp];
end

% --- Initial soil moisture and matric potential

for i = 1:soil.nsoi
   if (ityp == 0)
      soil.theta(i) = 0.10;
   elseif (ityp == 1)
      soil.theta(i) = 0.10;
   elseif (ityp == 2)
      soil.theta(i) = 0.24;
   end
   soil.psi(i) = matric_potential (soil.functions, params, soil.theta(i));
end

% --- Surface boundary condition: saturation (minus some small delta)

soil.theta0 = theta_sat - 1.0e-03;
if (ityp == 1)
   soil.theta0 = 0.267;
end
soil.psi0 = matric_potential (soil.functions, params, soil.theta0);

% --- Convergence criterion for delta_psi and water_balance_error

dpsi_tolerance = 1.0e-06;
water_balance_error = 1.0e-06;

% --- Time step (seconds)

dt = 10;
if (ityp == 1)
   dt = 5;
end

% --- Length of simulation (number of time steps)

% Hornberger & Wiberg: 15, 30, or 60 minutes
if (ityp == 0)
%  ntim = 15 * 60 / dt;
%  ntim = 30 * 60 / dt;
   ntim = 60 * 60 / dt;
end

% Haverkamp et al. (1977) - sand: duration is in hours
if (ityp == 1)
%  ntim = 0.05 * 3600 / dt;
%  ntim = 0.1 * 3600 / dt;
%  ntim = 0.2 * 3600 / dt;
%  ntim = 0.3 * 3600 / dt;
%  ntim = 0.4 * 3600 / dt;
   ntim = 0.8 * 3600 / dt;
end

% Haverkamp et al. (1977) - Yolo light clay: duration is in seconds
if (ityp == 2)
%  ntim = 1.0e4 / dt;
%  ntim = 1.0e5 / dt;
%  ntim = 5.0e5 / dt;
   ntim = 1.0e6 / dt;
end

% --- Initialize accumulators for water balance check

sum_in = 0;
sum_out = 0;
sum_store = 0;

% --- Time stepping loop: NTIM iterations with a time step of DT seconds

for itim = 1:ntim

   % Hour of day

   hour = itim * (dt/86400 * 24);
   fprintf('hour = %8.3f \n',hour)

   % Calculate soil moisture

   [soil] = picard (soil, params, dt, dpsi_tolerance, water_balance_error);

   % Sum fluxes for relative mass balance error

   sum_in = sum_in + abs(soil.Q0) * dt;
   sum_out = sum_out + abs(soil.QN) * dt;
   sum_store = sum_store + soil.dtheta;

   % Save cumulative infiltration

   xout(itim) = hour;
   yout(itim) = sum_in;

end

% --- Print mass balance

fprintf('infiltration (cm) = %8.3f \n',sum_in)
fprintf('drainage (cm) = %8.3f \n',sum_out)
fprintf('storage (cm) = %8.3f \n',sum_store)
relerr = ((sum_in - sum_out) - sum_store) / (sum_in - sum_out) * 100;
fprintf('mass balance error (percent) = %8.3f \n',relerr)

% --- Graph data

plot(soil.theta,soil.z)
xlabel('Volumetric moisture content')
ylabel('Depth (cm)')

% --- Write soil moisture as output file

A = [soil.theta; soil.z];
fileID = fopen('data1.txt','w');
fprintf(fileID,'%12s %12s\n','theta','z');
fprintf(fileID,'%12.3f %12.3f\n', A);
fclose(fileID);

% --- Write cumulative infiltration as output file

B = [xout; yout];
fileID = fopen('data2.txt','w');
fprintf(fileID,'%12s %12s\n','hour','infil');
fprintf(fileID,'%12.5f %12.5f\n', B);
fclose(fileID);
