% Supplemental program 7.1

% --- Physical constants

physcon.vkc = 0.4;                              % von Karman constant
physcon.grav = 9.80665;                         % Gravitational acceleration (m/s2)
physcon.tfrz = 273.15;                          % Freezing point of water (K)
physcon.sigma = 5.67e-08;                       % Stefan-Boltzmann constant (W/m2/K4)
physcon.mmdry = 28.97 / 1000;                   % Molecular mass of dry air (kg/mol)
physcon.mmh2o = 18.02 / 1000;                   % Molecular mass of water (kg/mol)
physcon.cpd = 1005.0;                           % Specific heat of dry air at constant pressure (J/kg/K)
physcon.cpw = 1846.0;                           % Specific heat of water vapor at constant pressure (J/kg/K)
physcon.rgas = 8.31446;                         % Universal gas constant (J/K/mol)
physcon.cwat = 4188.0;                          % Specific heat of water (J/kg/K)
physcon.cice = 2117.27;                         % Specific heat ice (J/kg/K)
physcon.rhowat = 1000.0;                        % Density of water (kg/m3)
physcon.rhoice = 917.0;                         % Density of ice (kg/m3)
physcon.cvwat = physcon.cwat * physcon.rhowat;  % Heat capacity of water (J/m3/K)
physcon.cvice = physcon.cice * physcon.rhoice;  % Heat capacity of ice (J/m3/K)
physcon.tkwat = 0.57;                           % Thermal conductivity of water (W/m/K)
physcon.tkice = 2.29;                           % Thermal conductivity of ice (W/m/K)
physcon.hfus = 0.3337e6;                        % Heat of fusion for water at 0 C (J/kg)
physcon.hvap = 2.501e6;                         % Latent heat of evaporation (J/kg)
physcon.hsub = physcon.hfus + physcon.hvap;     % Latent heat of sublimation (J/kg)

% --- Initialize soil texture variables

% Soil texture classes (Cosby et al. 1984. Water Resources Research 20:682-690)

%  1: sand
%  2: loamy sand
%  3: sandy loam
%  4: silty loam
%  5: loam
%  6: sandy clay loam
%  7  silty clay loam
%  8: clay loam
%  9: sandy clay
% 10: silty clay
% 11: clay

soilvar.silt = [ 5.0, 12.0, 32.0, 70.0, 39.0, 15.0, 56.0, 34.0,  6.0, 47.0, 20.0]; % Percent silt
soilvar.sand = [92.0, 82.0, 58.0, 17.0, 43.0, 58.0, 10.0, 32.0, 52.0,  6.0, 22.0]; % Percent sand
soilvar.clay = [ 3.0,  6.0, 10.0, 13.0, 18.0, 27.0, 34.0, 34.0, 42.0, 47.0, 58.0]; % Percent clay

% Volumetric soil water content at saturation (porosity)
% (Clapp and Hornberger. 1978. Water Resources Research 14:601-604)

soilvar.watsat = [0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.492, 0.482];

% --- Model run control parameters

dt = 1800;                          % Time step (seconds)
nday = 30;                          % Number of days to simulate, repeating the same diurnal cycle
soilvar.method = 'excess-heat';     % Phase change: use 'excess-heat' or 'apparent-heat-capacity'

  fluxvar.profiles = 'MOST';        % Use Monin-Obukhov similarity theory
% fluxvar.profiles = 'RSL';         % Use canopy coupling with roughness sublayer theory

  fluxvar.bucket = 'no_bucket';     % Soil wetness factor = 1
% fluxvar.bucket = 'use_bucket';    % Use bucket model hydrology for soil wetness factor

% --- Atmospheric forcing at a reference height

forcvar.zref = 30.0;      % Reference height (m)
tmean = 25.0;             % Mean daily air temperature (C)
trange = 10.0;            % Temperature range for diurnal cycle (C)
relhum = 70.0;            % Relative humidity (%)
forcvar.uref = 3.0;       % Wind speed at reference height (m/s)
forcvar.pref = 101325;    % Atmospheric pressure (Pa)
forcvar.rain = 0.0;       % Rainfall (kg H2O/m2/s)
forcvar.snow = 0.0;       % Snowfall (kg H2O/m2/s)

doy = 182.0;              % Day of year (1 to 365) for solar radiation
lat = 40.0 * pi/180;      % Latitude (degrees -> radians) for solar radiation
solcon = 1364.0;          % Solar constant (W/m2)

% --- Site characteristics

vis = 1; nir = 2;               % Waveband indices for visible and near-infrared
alb_surf(vis) = 0.10;           % Snow-free surface albedo for visible waveband (-)
alb_surf(nir) = 0.20;           % Snow-free surface albedo for near-infrared waveband (-)
alb_snow(vis) = 0.95;           % Snow albedo for visible waveband (-)
alb_snow(nir) = 0.70;           % Snow albedo for near-infrared waveband (-)
surfvar.emiss = 0.98;           % Surface emissivity (dimensionless)
snow_mask = 100.0;              % Snow albedo masking depth (kg H2O/m2)

soilvar.soil_texture = 5;       % Soil texture class

bucket.soil_water_max = 150.0;  % Maximum soil water (kg H2O/m2)
bucket.soil_beta_max = 0.75;    % Soil water at which soil_beta = 1 (fraction of soil_water_max)

surfvar.hc = 20.0;              % Canopy height (m)

%surfvar.hc = 0.5;
%forcvar.zref = 10.5;

surfvar.LAI = 5.0;              % Leaf area index (m2/m2)

% For Harman and Finnigan (2007, 2008) roughness sublayer (RSL)

lad = surfvar.LAI / surfvar.hc; % Leaf area density (m2/m3)
cd = 0.20;                      % Leaf drag coefficient (dimensionless)
surfvar.Lc = 1 / (cd * lad);    % Canopy density length scale (m)
surfvar.rc = 0.2;               % Leaf Nusselt number (heat) or Stanton number (scalar)

% For Monin-Obukhov parameterization (MOST)

fluxvar.disp = 0.67 * surfvar.hc;   % Displacement height (m)
fluxvar.z0m = 0.13 * surfvar.hc;    % Roughness length for momentum (m)
fluxvar.z0c = 0.10 * fluxvar.z0m;   % Roughness length for scalars (m)

% --- Soil variables

% Number of layers in soil profile

soilvar.nsoi = 10;

% Soil layer thickness (m)

soilvar.dz = [0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058];

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

% --- Initial conditions

% Initial soil temperature (K) and unfrozen and frozen water (kg H2O/m2)

for i = 1:soilvar.nsoi

   % Temperature

   soilvar.tsoi(i) = tmean + physcon.tfrz;

   % Soil water at saturation (kg H2O/m2)

   h2osoi_sat = soilvar.watsat(soilvar.soil_texture) * physcon.rhowat * soilvar.dz(i);

   % Actual water content is some fraction of saturation. These are only used for soil
   % thermal properties and phase change. Note the inconsistency with the use of soil
   % water in the bucket model to calculate the soil wetness factor.

   satfrac = 0.85;
   if (soilvar.tsoi(i) > physcon.tfrz)
      soilvar.h2osoi_ice(i) = 0;
      soilvar.h2osoi_liq(i) = satfrac * h2osoi_sat;
   else
      soilvar.h2osoi_liq(i) = 0;
      soilvar.h2osoi_ice(i) = satfrac * h2osoi_sat;
   end

end

% Initial surface temperature (K) and vapor pressure (Pa)

fluxvar.tsrf = soilvar.tsoi(1);
[esat, desat] = satvap (fluxvar.tsrf-physcon.tfrz);
fluxvar.esrf = esat;

% Bucket model snow and soil water (kg H2O/m2)

bucket.snow_water = 0;
bucket.soil_water = bucket.soil_water_max;

% --- Time stepping loop

% Main loop is NTIM iterations per day with a time step of DT seconds.
% This is repeated NDAY times to spinup the model from arbitrary initial
% soil temperatures.

ntim = round(86400/dt);
for j = 1:nday
   fprintf('day = %6.0f\n',j)
   for i = 1:ntim

      % Hour of day (0 to 24)

      hour = i * (dt/86400 * 24);

      % Air temperature (K): use a sine wave with max (tmean + 1/2 trange) at 1400
      % and min (tmean - 1/2 trange) at 0200

      forcvar.tref = tmean + 0.5 * trange * sin(2*pi/24 * (hour-8)) + physcon.tfrz;

      % Vapor pressure (Pa) using constant relative humidity

      [esat, desat] = satvap (forcvar.tref-physcon.tfrz);
      forcvar.eref = (relhum / 100) * esat;

      % Derived quantities
      % forcvar.thref   ! Potential temperature at reference height (K)
      % forcvar.qref    ! Specific humidity at reference height (kg/kg)
      % forcvar.thvref  ! Virtual potential temperature at reference height (K)
      % forcvar.rhomol  ! Molar density at reference height (mol/m3)
      % forcvar.rhoair  ! Air density at reference height (kg/m3)
      % forcvar.mmair   ! Molecular mass of air at reference height (kg/mol)
      % forcvar.cpair   ! Specific heat of air at constant pressure, at reference height (J/mol/K)

      forcvar.thref =  forcvar.tref + 0.0098 * forcvar.zref;
      forcvar.qref = physcon.mmh2o / physcon.mmdry * forcvar.eref / ...
         (forcvar.pref - (1 - physcon.mmh2o/physcon.mmdry) * forcvar.eref);
      forcvar.thvref = forcvar.thref * (1 + 0.61 * forcvar.qref);
      forcvar.rhomol = forcvar.pref / (physcon.rgas * forcvar.tref);
      forcvar.rhoair = forcvar.rhomol * physcon.mmdry * ...
         (1 - (1 - physcon.mmh2o/physcon.mmdry) * forcvar.eref / forcvar.pref);
      forcvar.mmair = forcvar.rhoair / forcvar.rhomol;
      forcvar.cpair = physcon.cpd * (1 + (physcon.cpw/physcon.cpd - 1) * forcvar.qref);     % J/kg/K
      forcvar.cpair = forcvar.cpair * forcvar.mmair;                                        % J/mol/K

      % Solar radiation (W/m2)
      % doy        ! Day of year (1 to 365)
      % lat        ! Latitude (radians)
      % decl       ! Declination (radians): Brock, T.D. (1981) Calculating solar radiation
      %            ! for ecological studies. Ecological Modelling 14:1-19
      % hour_angle ! Solar hour angle (radians)
      % coszen     ! Cosine of solar zenith angle
      % rv         ! Radius vector: Brock, T.D. 1981. Calculating solar radiation
      %            ! for ecological studies. Ecological Modelling 14:1-19
      % soltoa     ! Solar radiation on horizontal surface at top of atmosphere (W/m2)
      % tau_atm    ! Atmospheric transmission coefficient
      % oam        ! Optical air mass
      % soldir     ! Direct beam solar radiation on horizontal surface (W/m2)
      % soldif     ! Diffuse solar radiation on horizontal surface (W/m2)
      % solrad     ! Total solar radiation on horizontal surface (W/m2)

      % Solar radiation at top of the atmosphere

      decl = 23.45 * sin((284+doy)/365*2*pi) * pi/180;
      hour_angle = 15 * (hour-12) * pi/180;
      coszen = max(cos(lat)*cos(decl)*cos(hour_angle) + sin(lat)*sin(decl), 0);
      rv = 1 / sqrt(1 + 0.033*cos(doy/365*2*pi));
      soltoa = solcon / rv^2 * coszen;
      
      % Clear sky atmospheric attenuation: Gates, D.M. (1980) Biophysical Ecology, page 110, 115

      tau_atm = 0.5;
      oam = 1 / max(coszen, 0.04);
      soldir = soltoa * tau_atm^oam;                       % Clear sky direct beam
      soldif = soltoa * (0.271 - 0.294 * tau_atm^oam);     % Clear sky diffuse
      forcvar.solrad = soldir + soldif;                    % Total at surface

      % Longwave radiation (W/m2)

      forcvar.lwdown = (0.398e-05 * forcvar.tref^2.148) * physcon.sigma * forcvar.tref^4;

      % Effective surface albedo is weighted combination of snow-free and
      % snow albedos

      fsno = bucket.snow_water / (bucket.snow_water + snow_mask);
      alb_eff(vis) = alb_surf(vis) * (1 - fsno) + alb_snow(vis) * fsno;
      alb_eff(nir) = alb_surf(nir) * (1 - fsno) + alb_snow(nir) * fsno;

      % Radiative forcing: absorbed solar + incident longwave. This partitions
      % solar radiation into 50% visible and 50% near-infrared wavebands.

      fluxvar.qa = (1-alb_eff(vis)) * 0.5*forcvar.solrad ...
      + (1-alb_eff(nir)) * 0.5*forcvar.solrad + surfvar.emiss * forcvar.lwdown;

      % Canopy conductance (mol/m2/s) - use a weighted average of sunlit and shaded leaves

      gcan_min = 0.05;                           % Minimum conductance (mol/m2/s)
      gcan_max = 0.2;                            % Maximum conductance (mol/m2/s)

      ext = 0.5 / max(coszen, 0.0001);           % Light extinction coefficient
      fsun = (1 - exp(-ext*surfvar.LAI)) / (ext*surfvar.LAI);  % Sunlit fraction of canopy
      if (soldir+soldif > 0)
         surfvar.gcan = (fsun * gcan_max + (1 - fsun) * gcan_min) * surfvar.LAI;
      else
         surfvar.gcan = gcan_min * surfvar.LAI;
      end

      % Thermal conductivity and heat capacity

     [soilvar] = soil_thermal_properties (physcon, soilvar);

      % Calculate the soil temperatures and surface fluxes

      [fluxvar, soilvar, bucket] = surface_fluxes (physcon, forcvar, surfvar, soilvar, fluxvar, bucket, dt);

      % Rainfall to equal evaporative loss (kg H2O/m2/s)

%     forcvar.rain = fluxvar.etflx * physcon.mmh2o;
      forcvar.rain = 0;

      % Bucket model hydrology

      [bucket] = bucket_hydrology (physcon, forcvar, fluxvar, bucket, dt);

      % Save data for graphics

      xhour(i) = hour;
      ytsrf(i) = fluxvar.tsrf - physcon.tfrz;
      ytref(i) = forcvar.tref - physcon.tfrz;
      yrnet(i) = fluxvar.rnet;
      yshflx(i) = fluxvar.shflx;
      ylhflx(i) = fluxvar.lhflx;
      ygsoi(i) = fluxvar.gsoi;
      yustar(i) = fluxvar.ustar;
      ygac(i) = fluxvar.gac * 100;
      ygcan(i) = surfvar.gcan * 100;

   end
end

% --- Write output files

A = [xhour; ytref; ytsrf; yrnet; yshflx; ylhflx; ygsoi; yustar; ygac; ygcan];
fileID = fopen('flux.txt','w');
fprintf(fileID,'%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n','hour','Ta','Ts','Rn','H','LE','G','ustar','gac','gcan');
fprintf(fileID,'%12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n', A);
fclose(fileID);

B = [soilvar.z; soilvar.tsoi];
fileID = fopen('tsoi.txt','w');
fprintf(fileID,'%12s %12s\n','depth','tsoi');
fprintf(fileID,'%12.3f %12.3f\n', B);
fclose(fileID);

% --- Make graph

plot(xhour,yrnet,'g-',xhour,yshflx,'r-',xhour,ylhflx,'b-',xhour,ygsoi,'r--',xhour,ygac,'m-')
axis([0 24 -100 600])
set(gca,'xTick',0:3:24)
set(gca,'yTick',-100:100:800)
title('Diurnal cycle')
xlabel('Time of day (hours)')
ylabel('Flux (W m^{-2})')
legend('R_n','H','\lambdaE','G','g_{ac}*100','Location','northwest')
