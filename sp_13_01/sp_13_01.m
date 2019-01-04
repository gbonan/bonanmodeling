% Supplemental program 13.1

% -------------------------------------------------------------------------
% Calculate leaf gas exchange coupled with the leaf energy budget for C3
% and C4 plants. Leaf temperature is calculated from the energy balance.
% Stomatal conductance is calculated from water-use efficiency optimization
% within the constraints imposed by plant hydraulics.
% -------------------------------------------------------------------------

% --- Waveband indices for visible and near-infrared

params.vis = 1; params.nir = 2;

% --- Physical constants

physcon.grav = 9.80665;               % Gravitational acceleration (m/s2)
physcon.tfrz = 273.15;                % Freezing point of water (K)
physcon.sigma = 5.67e-08;             % Stefan-Boltzmann constant (W/m2/K4)
physcon.mmdry = 28.97 / 1000;         % Molecular mass of dry air (kg/mol)
physcon.mmh2o = 18.02 / 1000;         % Molecular mass of water (kg/mol)
physcon.cpd = 1005;                   % Specific heat of dry air at constant pressure (J/kg/K)
physcon.cpw = 1846;                   % Specific heat of water vapor at constant pressure (J/kg/K)
physcon.rgas = 8.31446;               % Universal gas constant (J/K/mol)
physcon.visc0 = 13.3e-06;             % Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
physcon.Dh0 = 18.9e-06;               % Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
physcon.Dv0 = 21.8e-06;               % Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
physcon.Dc0 = 13.8e-06;               % Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)
physcon.denh2o = 1000;                % Density of liquid water (kg/m3)

% --- Set leaf physiology variables

% Photosynthetic pathway: 1 = C3. 0 = C4

leaf.c3psn = 1;

% Photosynthesis co-limitation: 0 = no. 1 = yes

leaf.colim = 1;

% Leaf physiological parameters

[leaf] = LeafPhysiologyParams (params, physcon, leaf);

% --- Set root parameters

[rootvar] = RootParams();

% --- Set site variables

% Leaf height (m)

flux.height = 15;

% Fine root biomass (g biomass / m2)

rootvar.biomass = 500;

% Canopy leaf area index (m2/m2)

flux.lai = 5;

% Soil texture classes
%  1: sand
%  2: loamy sand
%  3: sandy loam
%  4: silt loam
%  5: loam
%  6: sandy clay loam
%  7  silty clay loam
%  8: clay loam
%  9: sandy clay
% 10: silty clay
% 11: clay

soil.texture = 5;

% --- Set soil hydraulic parameters, soil depth, and rooting fraction

[soil] = SoilParams (soil);

% --- Set soil moisture

for j = 1:soil.nlevsoi
   soil.h2osoi_vol(j) = 0.5 * soil.watsat(j);
   soil.psi(j) = soil.psisat(j) * (soil.h2osoi_vol(j) / soil.watsat(j))^-soil.bsw(j);
end

% --- Hydraulic conductances

% Aboveground plant stem resistance for xylem-to-leaf (MPa.s.m2/mmol H2O)

flux.rplant = 1 / leaf.gplant;

% Calculate soil hydraulic resistance, weighted soil water potential and
% water uptake from each soil layer

[flux] = SoilResistance (physcon, leaf, rootvar, soil, flux);

% Leaf-specific conductance for soil-to-leaf (mmol H2O/m2 laef/s/MPa)

flux.lsc = 1 / (flux.rsoil + flux.rplant);

fprintf('rplant = %15.5f\n',flux.rplant)
fprintf('rsoil = %15.5f\n',flux.rsoil)
fprintf('lsc = %15.5f\n',flux.lsc)

% --- Atmospheric forcing

% Process sunlit or shaded leaf

  leaftype = 'sun';
% leaftype = 'shade';

% Atmospheric CO2 (umol/mol) and O2 (mmol/mol)

atmos.co2air = 380;
atmos.o2air = 0.209 * 1000;

% Air temperature (K) and relative humidity (%)

atmos.tair = physcon.tfrz + 30;
atmos.relhum = 60;

% Wind (m/s)
% u = 0.01_r8  ! Still air
% u = 0.1_r8   ! Calm - smoke rises vertically
% u = 1.0_r8   ! Light air - smoke drift indicates wind direction
% u = 2.5_r8   ! Light breeze - wind felt on skin, leaves rustle
% u = 5.0_r8   ! Gentle breeze - leaves constantly moving and light flag extended
% u = 10.0_r8  ! Moderate breeze

atmos.wind = 1;

% Atmospheric pressure (Pa)

atmos.patm = 101325;

% Vapor pressure (Pa) and specific humidity (kg/kg)

[esat, desat] = satvap ((atmos.tair-physcon.tfrz));
atmos.eair = esat * (atmos.relhum / 100);
atmos.qair = physcon.mmh2o / physcon.mmdry * atmos.eair / (atmos.patm - (1 - physcon.mmh2o/physcon.mmdry) * atmos.eair);

% Molar density (mol/m3)

atmos.rhomol = atmos.patm / (physcon.rgas * atmos.tair);

% Air density (kg/m3)

atmos.rhoair = atmos.rhomol * physcon.mmdry * (1 - (1 - physcon.mmh2o/physcon.mmdry) * atmos.eair / atmos.patm);

% Molecular mass of air (kg/mol)

atmos.mmair = atmos.rhoair / atmos.rhomol;

% Specific heat of air at constant pressure (J/mol/K)

atmos.cpair = physcon.cpd * (1 + (physcon.cpw/physcon.cpd - 1) * atmos.qair) * atmos.mmair;

% Atmospheric longwave radiation (W/m2)

atmos.irsky = 400;

% Solar radiation (W/m2)

switch leaftype
   case 'sun'
   fsds = 800;   % Sun leaf
   case 'shade'
   fsds = 300;   % Shade leaf
end
atmos.swsky(params.vis) = 0.5 * fsds;
atmos.swsky(params.nir) = 0.5 * fsds;

% --- Ground variables

ground.albsoi(params.vis) = 0.1;      % Soil albedo (visible waveband)
ground.albsoi(params.nir) = 0.2;      % Soil albedo (near-infrared waveband)
tg = atmos.tair;
ground.irgrd = physcon.sigma * tg^4;

% --- Radiation absorbed by leaf

% Solar radiation incident on leaf

flux.swinc(params.vis) = atmos.swsky(params.vis) * (1 + ground.albsoi(params.vis));
flux.swinc(params.nir) = atmos.swsky(params.nir) * (1 + ground.albsoi(params.nir));

% Solar radiation absorbed by leaf

flux.swflx(params.vis) = flux.swinc(params.vis) * (1 - leaf.rho(params.vis) - leaf.tau(params.vis));
flux.swflx(params.nir) = flux.swinc(params.nir) * (1 - leaf.rho(params.nir) - leaf.tau(params.nir));
flux.apar = flux.swflx(params.vis) * 4.6;

% Radiative forcing for leaf temperature calculation

flux.qa = flux.swflx(params.vis) + flux.swflx(params.nir) + leaf.emiss * (atmos.irsky + ground.irgrd);

% --- Leaf temperature, energy fluxes, photosynthesis, and stomatal conductance

% Initial estimate for leaf temperature

flux.tleaf = atmos.tair;

% Leaf water potential at beginning of time step (MPa)

flux.psi_leaf = -0.1;

% Time step (s) for change in leaf water potential

flux.dt = 30 * 60;

% Flux calculations

[flux] = LeafFluxes (physcon, atmos, leaf, flux);

% Print flux output

fprintf('dleaf = %15.5f\n',leaf.dleaf*100)          % m -> cm
fprintf('apar = %15.5f\n',flux.apar)
fprintf('tleaf = %15.5f\n',flux.tleaf-physcon.tfrz) % K -> oC
fprintf('qa = %15.5f\n',flux.qa)
fprintf('lhflx = %15.5f\n',flux.lhflx)
fprintf('etflx = %15.5f\n',flux.etflx*1000)         % mol H2O/m2/s -> mmol H2O/m2/s
fprintf('an = %15.5f\n',flux.an)
fprintf('gbh = %15.5f\n',flux.gbh)
fprintf('gs = %15.5f\n',flux.gs)
fprintf('psi_leaf = %15.5f\n',flux.psi_leaf)
