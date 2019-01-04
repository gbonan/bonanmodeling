% Supplemental program 5.3

% -------------------------------------------------------------------------
% Use implicit formulation with "excess heat" or "apparent heat capacity"
% to solve for soil temperatures with phase change in comparison with
% Neumann's analytical solution.
% -------------------------------------------------------------------------

% --- Physical constants in physcon structure

physcon.tfrz = 273.15;             % Freezing point of water (K)
physcon.rhowat = 1000;             % Density of water (kg/m3)
physcon.rhoice = 917;              % Density of ice (kg/m3)
physcon.hfus = 0.3337e6;           % Heat of fusion for water at 0 C (J/kg)

% --- Model run control parameters

dt = 3600;                         % Time step (seconds)
nday = 60;                         % Number of days
%soilvar.method = 'excess-heat';            % Use excess heat for phase change
soilvar.method = 'apparent-heat-capacity'; % Use apparent heat capacity for phase change

% --- Initialize soil layer variables

% Number of layers in soil profile

soilvar.nsoi = 60;

% Soil layer thickness (m)

for i = 1:soilvar.nsoi
   soilvar.dz(i) = 0.10;
end

% Soil depth (m) at i+1/2 interface between layers i and i+1 (negative distance from surface)

soilvar.z_plus_onehalf(1) = -soilvar.dz(1);
for i = 2:soilvar.nsoi
   soilvar.z_plus_onehalf(i) = soilvar.z_plus_onehalf(i-1) - soilvar.dz(i);
end

% Soil depth (m) at center of layer i (negative distance from surface)

soilvar.z(1) = 0.5 * soilvar.z_plus_onehalf(1);
for i = 2:soilvar.nsoi
   soilvar.z(i) = 0.5 * (soilvar.z_plus_onehalf(i-1) + soilvar.z_plus_onehalf(i));
end

% Thickness between between z(i) and z(i+1)

for i = 1:soilvar.nsoi-1
   soilvar.dz_plus_onehalf(i) = soilvar.z(i) - soilvar.z(i+1);
end
soilvar.dz_plus_onehalf(soilvar.nsoi) = 0.5 * soilvar.dz(soilvar.nsoi);

% Initial soil temperature (K)

for i = 1:soilvar.nsoi
   soilvar.tsoi(i) = physcon.tfrz + 2;
end

% Initial unfrozen and frozen water (kg H2O/m2)

for i = 1:soilvar.nsoi
   if (soilvar.tsoi(i) > physcon.tfrz)
      soilvar.h2osoi_ice(i) = 0;
      soilvar.h2osoi_liq(i) = 0.187 * 1770 * soilvar.dz(i);
   else
      soilvar.h2osoi_liq(i) = 0;
      soilvar.h2osoi_ice(i) = 0.187 * 1770 * soilvar.dz(i);
   end
end

% Soil layers for output

k1 = 3;
k2 = 6;
k3 = 9;
k4 = 12;

% --- Time stepping loop to increment soil temperature

% Counter for output file

m = 1;

% Save initial data for output files

iday_out(m) = 0;
d0c_out(m) = 0;
z1_out(m) = -soilvar.z(k1) * 100;
tsoi1_out(m) = soilvar.tsoi(k1) - physcon.tfrz;
z2_out(m) = -soilvar.z(k2) * 100;
tsoi2_out(m) = soilvar.tsoi(k2) - physcon.tfrz;
z3_out(m) = -soilvar.z(k3) * 100;
tsoi3_out(m) = soilvar.tsoi(k3) - physcon.tfrz;
z4_out(m) = -soilvar.z(k4) * 100;
tsoi4_out(m) = soilvar.tsoi(k4) - physcon.tfrz;

% Main loop is NTIM iterations per day with a time step of DT seconds.
% This is repeated NDAY times.

ntim = round(86400/dt);
for iday = 1:nday
   for itim = 1:ntim

      tsurf = -10 + physcon.tfrz;

      % Thermal conductivity and heat capacity

      [soilvar] = soil_thermal_properties (physcon, soilvar);

      % Soil temperatures

      [soilvar] = soil_temperature (physcon, soilvar, tsurf, dt);

      % Calculate depth to freezing isotherm (0 C)

      d0c = 0; % Depth of 0 C isotherm (m)

      switch soilvar.method
         case 'excess-heat'
         num_z = 0;
         sum_z = 0;

         % Average depth of soil layers with TSOI = TFRZ

         for i = 1:soilvar.nsoi
            if (abs(soilvar.tsoi(i)-physcon.tfrz) < 1.e-03)
               % Number of soil layers for averaging
               num_z = num_z + 1;
               % Sum of soil layer depths for averaging
               sum_z = sum_z + soilvar.z(i);
            end
         end

         if (num_z > 0)
            d0c = sum_z / num_z;
         end

         case 'apparent-heat-capacity'
         % Linear interpolation between soil layers

         for i = 2:soilvar.nsoi
            if (soilvar.tsoi(i-1) <= physcon.tfrz & soilvar.tsoi(i) > physcon.tfrz)
               % slope for linear interpolation
               b = (soilvar.tsoi(i) - soilvar.tsoi(i-1)) / (soilvar.z(i) - soilvar.z(i-1));
               % intercept
               a = soilvar.tsoi(i) - b * soilvar.z(i);
               % depth of 0C isotherm
               d0c = (physcon.tfrz - a) / b;
            end
         end
      end

      % Save data for output files

      if (itim == ntim)
         m = m + 1;
         iday_out(m) = iday;
         d0c_out(m) = -d0c * 100;
         z1_out(m) = -soilvar.z(k1) * 100;
         tsoi1_out(m) = soilvar.tsoi(k1) - physcon.tfrz;
         z2_out(m) = -soilvar.z(k2) * 100;
         tsoi2_out(m) = soilvar.tsoi(k2) - physcon.tfrz;
         z3_out(m) = -soilvar.z(k3) * 100;
         tsoi3_out(m) = soilvar.tsoi(k3) - physcon.tfrz;
         z4_out(m) = -soilvar.z(k4) * 100;
         tsoi4_out(m) = soilvar.tsoi(k4) - physcon.tfrz;
      end

   end
end

A = [iday_out; d0c_out; z1_out; tsoi1_out; z2_out; tsoi2_out; z3_out; tsoi3_out; z4_out; tsoi4_out];
fileID = fopen('data_numerical.txt','w');
fprintf(fileID,'%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n','day','z0C','z1','t1','z2','t2','z3','t3','z4','t4');
fprintf(fileID,'%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n', A);
fclose(fileID);

% Analytical solution for Neumann problem (Lunardini 1981)

[dummy] = neumann;
