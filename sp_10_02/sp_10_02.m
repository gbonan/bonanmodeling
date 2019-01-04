% Supplemental program 10.2

% -------------------------------------------------------------------------
% Calculate leaf temperature for different radiative forcing, wind speed,
% and stomatal conductance
% -------------------------------------------------------------------------

% Pass variables to other routines as global variables

global tfrz sigma g visc0 Dh0 Dv0 Dc0 mmh2o
global tair pref eair wind cpair rhomol
global gsw gbw gbh gbc emleaf dleaf
global tleaf qa rn lwrad sh lh

% Physical constants

tfrz = 273.15;            % Freezing point of water (K)
sigma = 5.67e-08;         % Stefan-Boltzmann constant (W/m2/K4)
g = 9.80616;              % Gravitational acceleration (m/s2)
visc0 = 13.3e-06;         % Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
Dh0 = 18.9e-06;           % Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
Dv0 = 21.8e-06;           % Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
Dc0 = 13.8e-06;           % Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)
rgas = 8.31447;           % Universal gas constant (J/K/mol)
mmdry = 28.966 / 1000;    % Molecular mass of dry air (kg/mol)
mmh2o = 18.016 / 1000;    % Molecular mass of water vapor (kg/mol)
cpd = 1004.64;            % Specific heat of dry air at constant pressure (J/kg/K)
cpw = 1810;               % Specific heat of water vapor at constant pressure (J/kg/K)

% Waveband indices: 1 = visible. 2 = near-infrared

vis = 1;
nir = 2;

% Number of data points to simulate: 4 wind speeds * 2 stomatal conductance * 2 solar radiation

num = 4 * 2 * 2;

% Input variables

for p = 1:num

   tair(p) = tfrz + 35;                  % Air temperature (K)
   relhum(p) = 50;                       % Relative humidity (%)
   pref(p) = 101325;                     % Air pressure (Pa)
   irsky(p) = 300;                       % Atmospheric longwave radiation (W/m2)
   irgrd(p) = sigma * (tfrz + 39.95)^4;  % Ground longwave radiation (W/m2)

   % Set wind speed (m/s) to specific values

   if (p == 1 | p == 5 | p == 9 | p == 13)
      wind(p) = 0.01;  % Still air
   end
   if (p == 2 | p == 6 | p == 10 | p == 14)
      wind(p) = 0.1;   % Calm - smoke rises vertically
   end
   if (p == 3 | p == 7 | p == 11 | p == 15)
      wind(p) = 1.0;   % Light air - smoke drift indicates wind direction
   end
   if (p == 4 | p == 8 | p == 12 | p == 16)
      wind(p) = 5.0;   % Gentle breeze - leaves constantly moving and light flag extended
   end

   % Solar radiation (W/m2) for visible and near-infrared wavebands

   if (p <= 8)
      swsky(p,vis) = 0.5 * 1100;         % Full sun
      swsky(p,nir) = 0.5 * 1100;
   else
      swsky(p,vis) = 0.5 *  550;         % Cloudy
      swsky(p,nir) = 0.5 *  550;
   end

   % Albedo of ground surface for visible and near-infrared wavebands

   albsoi(p,vis) = 0.1;
   albsoi(p,nir) = 0.2;

   % Leaf input

   rhol(p,vis) = 0.1;   % Leaf reflectance
   rhol(p,nir) = 0.4;
   taul(p,vis) = 0.1;   % Leaf transmittance
   taul(p,nir) = 0.4;

   emleaf(p) = 0.98;     % Leaf emissivity
   dleaf(p) = 0.05;      % Leaf dimension (m)

   % Leaf stomatal conductance (mol H2O/m2/s)

   if ((p >= 1 & p <=4) | (p >= 9 & p <= 12))
      gsw(p) = 0;            % Only longwave and sensible heat. No latent heat
   end
   if ((p >= 5 & p <=8 | p >= 13 & p <= 16))
      gsw(p) = 0.4;          % Longwave, sensible heat, and latent heat
   end

end

% Derived quantities

for p = 1:num

   % esat    ! Saturation vapor pressure of air (Pa)
   % eair    ! Vapor pressure of air (Pa)
   % qair    ! Specific humidity (kg/kg)
   % rhomol  ! Molar density (mol/m3)
   % rhoair  ! Air density (kg/m3)
   % mmair   ! Molecular mass of air (kg/mol)
   % cpair   ! Specific heat of air at constant pressure (J/mol/K)

   [esat, desat] = satvap (tair(p)-tfrz);
   eair(p) = (relhum(p) / 100) * esat;
   qair(p) = mmh2o / mmdry * eair(p) / (pref(p) - (1 - mmh2o/mmdry) * eair(p));
   rhomol(p) = pref(p) / (rgas * tair(p));
   rhoair(p) = rhomol(p) * mmdry * (1 - (1 - mmh2o/mmdry) * eair(p) / pref(p));
   mmair(p) = rhoair(p) / rhomol(p);
   cpair(p) = cpd * (1 + (cpw/cpd - 1) * qair(p)) * mmair(p);

end

% Radiative forcing (W/m2): absorbed solar + absorbed longwave radiation

for p = 1:num

   swinc(p,vis) = swsky(p,vis) * (1 + albsoi(p,vis));
   swinc(p,nir) = swsky(p,nir) * (1 + albsoi(p,nir));

   qa(p) = swinc(p,vis) * (1 - rhol(p,vis) - taul(p,vis)) ...
         + swinc(p,nir) * (1 - rhol(p,nir) - taul(p,nir)) + emleaf(p) * (irsky(p) + irgrd(p));

end

% Leaf temperature and fluxes

for p = 1:num

   rn(p) = 0;           % Leaf net radiation (W/m2)
   lwrad(p) = 0;        % Longwave radiation emitted from leaf (W/m2)
   sh(p) = 0;           % Leaf sensible heat flux (W/m2)
   lh(p) = 0;           % Leaf latent heat flux (W/m2)
   gbh(p) = 0;          % Leaf boundary layer conductance, heat (mol/m2/s)
   gbw(p) = 0;          % Leaf boundary layer conductance, H2O (mol H2O/m2/s)
   gbc(p) = 0;          % Leaf boundary layer conductance, CO2 (mol CO2/m2/s)
   tleaf(p) = tair(p);  % Initial estimate leaf temperature (K)

   % Solve for leaf temperature and fluxes. Need to iterate because boundary
   % layer conductances depend on tleaf (for free convection)

   niter = 0;
   delta = 1e36;

   while (niter <= 100 & abs(delta) > 1e-06)

      % Increment iteration counter

      niter = niter + 1;

      % Save temperature from previous iteration

      tleaf_old = tleaf(p);

      % Leaf boundary layer conductances

      [x] = leaf_boundary_layer (p);

      % Leaf temperature and energy fluxes

      [x] = leaf_temperature (p);

      % Change in leaf temperature

      delta = tleaf(p) - tleaf_old;

   end
end

tleaf = tleaf - tfrz;
fprintf('qa = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',qa(1),qa(2),qa(3),qa(4),qa(5),qa(6),qa(7),qa(8))
fprintf('gsw = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',gsw(1),gsw(2),gsw(3),gsw(4),gsw(5),gsw(6),gsw(7),gsw(8))
fprintf('wind = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',wind(1),wind(2),wind(3),wind(4),wind(5),wind(6),wind(7),wind(8))
fprintf('tleaf = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',tleaf(1),tleaf(2),tleaf(3),tleaf(4),tleaf(5),tleaf(6),tleaf(7),tleaf(8))
fprintf('lwrad = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',lwrad(1),lwrad(2),lwrad(3),lwrad(4),lwrad(5),lwrad(6),lwrad(7),lwrad(8))
fprintf('sh = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',sh(1),sh(2),sh(3),sh(4),sh(5),sh(6),sh(7),sh(8))
fprintf('lh = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',lh(1),lh(2),lh(3),lh(4),lh(5),lh(6),lh(7),lh(8))

fprintf(' \n')

fprintf('qa = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',qa(9),qa(10),qa(11),qa(12),qa(13),qa(14),qa(15),qa(16))
fprintf('gsw = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',gsw(9),gsw(10),gsw(11),gsw(12),gsw(13),gsw(14),gsw(15),gsw(16))
fprintf('wind = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',wind(9),wind(10),wind(11),wind(12),wind(13),wind(14),wind(15),wind(16))
fprintf('tleaf = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',tleaf(9),tleaf(10),tleaf(11),tleaf(12),tleaf(13),tleaf(14),tleaf(15),tleaf(16))
fprintf('lwrad = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',lwrad(9),lwrad(10),lwrad(11),lwrad(12),lwrad(13),lwrad(14),lwrad(15),lwrad(16))
fprintf('sh = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',sh(9),sh(10),sh(11),sh(12),sh(13),sh(14),sh(15),sh(16))
fprintf('lh = %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',lh(9),lh(10),lh(11),lh(12),lh(13),lh(14),lh(15),lh(16))
