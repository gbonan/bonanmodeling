% Supplemental program 10.1

% -----------------------------------------
% Calculate leaf boundary layer conductance
% -----------------------------------------

% Physical constants

tfrz = 273.15;                       % Freezing point of water (K)
g = 9.80616;                         % Gravitational acceleration (m/s2)
rgas = 8.31447;                      % Universal gas constant (J/K/mol)
visc0 = 13.3e-06;                    % Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
Dh0 = 18.9e-06;                      % Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
Dv0 = 21.8e-06;                      % Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
Dc0 = 13.8e-06;                      % Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)

% Calculations are for 100 data points that vary in wind speed and/or leaf temperature

num = 100;

% Input variables

for p = 1:num
   tair(p) = tfrz + 15;                      % Air temperature (K)
   pref(p) = 101325;                         % Pressure (Pa)
   rhomol(p) = pref(p) / (rgas * tair(p));   % Molar density (mol/m3)
   wind(p) = p / 10;                         % Wind speed (m/s): 0.1 -> 10 m/s by 0.1
   dleaf(p) = 0.05;                          % Leaf dimension (m)
   tleaf(p) = tair(p) + p / 10;              % tleaf - tair: 0.1 -> 10 K by 0.1 (for free convection)
end

for p = 1:num

   % Adjust diffusivity for temperature and pressure

   fac = 101325 / pref(p) * (tair(p) / tfrz)^1.81;

   visc = visc0 * fac;    % Kinematic viscosity (m2/s)
   Dh = Dh0 * fac;        % Molecular diffusivity, heat (m2/s)
   Dv = Dv0 * fac;        % Molecular diffusivity, H2O (m2/s)
   Dc = Dc0 * fac;        % Molecular diffusivity, CO2 (m2/s)

   % Dimensionless numbers

   Re(p) = wind(p) * dleaf(p) / visc;   % Reynolds number
   Pr  = visc / Dh;                     % Prandtl number
   Scv = visc / Dv;                     % Schmidt number for H2O
   Scc = visc / Dc;                     % Schmidt number for CO2
   Gr = g * dleaf(p)^3 * max(tleaf(p) - tair(p), 0) / (tair(p) * visc * visc); % Grashof number

   % Forced convection - laminar flow

   b1 = 1.5;                                     % Empirical correction factor for Nu and Sh

   Nu  = b1 * 0.66 *  Pr^0.33 * Re(p)^0.5;       % Nusselt number
   Shv = b1 * 0.66 * Scv^0.33 * Re(p)^0.5;       % Sherwood number, H2O
   Shc = b1 * 0.66 * Scc^0.33 * Re(p)^0.5;       % Sherwood number, CO2

   gbh1(p) = Dh * Nu / dleaf(p) * rhomol(p);     % Boundary layer conductance, heat (mol/m2/s)
   gbw1(p) = Dv * Shv / dleaf(p) * rhomol(p);    % Boundary layer conductance, H2O (mol/m2/s)
   gbc1(p) = Dc * Shc / dleaf(p) * rhomol(p);    % Boundary layer conductance, CO2 (mol/m2/s)

   % Forced convection - turbulent flow

   Nu  = b1 * 0.036 *  Pr^0.33 * Re(p)^0.8;      % Nusselt number
   Shv = b1 * 0.036 * Scv^0.33 * Re(p)^0.8;      % Sherwood number, H2O
   Shc = b1 * 0.036 * Scc^0.33 * Re(p)^0.8;      % Sherwood number, CO2

   gbh2(p) = Dh * Nu / dleaf(p) * rhomol(p);     % Boundary layer conductance, heat (mol/m2/s)
   gbw2(p) = Dv * Shv / dleaf(p) * rhomol(p);    % Boundary layer conductance, H2O (mol/m2/s)
   gbc2(p) = Dc * Shc / dleaf(p) * rhomol(p);    % Boundary layer conductance, CO2 (mol/m2/s)

   % Free convection

   Nu  = 0.54 *  Pr^0.25 * Gr^0.25;              % Nusselt number
   Shv = 0.54 * Scv^0.25 * Gr^0.25;              % Sherwood number, H2O
   Shc = 0.54 * Scc^0.25 * Gr^0.25;              % Sherwood number, CO2

   gbh3(p) = Dh * Nu / dleaf(p) * rhomol(p);     % Boundary layer conductance, heat (mol/m2/s)
   gbw3(p) = Dv * Shv / dleaf(p) * rhomol(p);    % Boundary layer conductance, H2O (mol/m2/s)
   gbc3(p) = Dc * Shc / dleaf(p) * rhomol(p);    % Boundary layer conductance, CO2 (mol/m2/s)

end

% Make graph

plot(wind,gbh1,'b-',wind,gbh2,'r-')
title('Heat')
xlabel('Wind speed (m s^{-1})')
ylabel('Boundary layer conductance (mol m^{-2} s^{-1})')
legend('Laminar','Turbulent','Location','best')

% Write output to files

A = [wind; Re; gbh1];
fileID = fopen('forced-laminar.dat','w');
fprintf(fileID,'%10s %10s %10s\n','wind','Re','gbh');
fprintf(fileID,'%10.1f %10.2f %10.4f\n', A);
fclose(fileID);

B = [wind; Re; gbh2];
fileID = fopen('forced-turbulent.dat','w');
fprintf(fileID,'%10s %10s %10s\n','wind','Re','gbh');
fprintf(fileID,'%10.1f %10.2f %10.4f\n', B);
fclose(fileID);

C = [tleaf-tair; gbh3];
fileID = fopen('free-convection.dat','w');
fprintf(fileID,'%10s %10s\n','Tl-Ta','gbh');
fprintf(fileID,'%10.1f %10.4f\n', C);
fclose(fileID);
