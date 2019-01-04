function [fluxvar, soilvar, bucket] = surface_fluxes (physcon, forcvar, surfvar, soilvar, fluxvar, bucket, dt)

% Calculate soil temperatures and surface fluxes. This routine uses the
% current estimate of surface temperature and vapor pressure to solve
% for the Obukhov length. This then provides the aerodynamic conductance,
% which is used in the surface temperature and flux calculation.

% ------------------------------------------------------
% Input
%   dt                  ! Time step (s)
%   physcon.tfrz        ! Freezing point of water (K)
%   physcon.mmh2o       ! Molecular mass of water (kg/mol)
%   physcon.hvap        ! Latent heat of evaporation (J/kg)
%   physcon.hsub        ! Latent heat of sublimation (J/kg)
%   physcon.sigma       ! Stefan-Boltzmann constant (W/m2/K4)
%   forcvar.thref       ! Potential temperature at reference height (K)
%   forcvar.uref        ! Wind speed at reference height (m/s)
%   forcvar.eref        ! Vapor pressure at reference height (Pa)
%   forcvar.pref        ! Atmospheric pressure (Pa)
%   forcvar.cpair       ! Specific heat of air at constant pressure, at reference height (J/mol/K)
%   forcvar.rhomol      ! Molar density at reference height (mol/m3)
%   surfvar.emiss       ! Surface emissivity
%   surfvar.gcan        ! Canopy conductance (mol/m2/s)
%   soilvar.method      ! Use excess heat or apparent heat capacity for phase change
%   soilvar.nsoi        ! Number of soil layers
%   soilvar.dz          ! Soil layer thickness (m)
%   soilvar.cv          ! Heat capacity (J/m3/K)
%   fluxvar.profiles    ! Use MOST or RSL for flux-profiles
%   fluxvar.qa          ! Radiative forcing (W/m2)
%   fluxvar.bucket      ! Use bucket model hydrology soil wetness factor
%   bucket.snow_water   ! Snow water (kg H2O/m2)
%   bucket.soil_water   ! Soil water (kg H2O/m2)
%   bucket.soil_water_max ! Maximum soil water (kg H2O/m2)
%   bucket.soil_beta_max  ! Soil water at which soil_beta = 1 (fraction of soil_water_max)
%
% Input/output
%   fluxvar.tsrf        ! Surface temperature (K)
%   fluxvar.esrf        ! Surface vapor pressure (Pa)
%   soilvar.tsoi        ! Soil temperature (K)
%   soilvar.h2osoi_liq  ! Unfrozen water, liquid (kg H2O/m2)
%   soilvar.h2osoi_ice  ! Frozen water, ice (kg H2O/m2)
%
% Output
%   fluxvar.rnet        ! Net radiation (W/m2)
%   fluxvar.lwrad       ! Emitted longwave radiation (W/m2)
%   fluxvar.shflx       ! Sensible heat flux (W/m2)
%   fluxvar.lhflx       ! Latent heat flux (W/m2)
%   fluxvar.etflx       ! Evapotranspiration (mol H2O/m2/s)
%   fluxvar.gsoi        ! Soil energy flux (W/m2)
%   fluxvar.gsno        ! Snow melt energy flux (W/m2)
%   fluxvar.gam         ! Aerodynamic conductance for momentum (mol/m2/s)
%   fluxvar.gac         ! Aerodynamic conductance for scalars (mol/m2/s)
%   soilvar.hfsoi       ! Soil phase change energy flux (W/m2)
%   bucket.soil_beta    ! Soil wetness factor for evapotranspiration (-)
%   bucket.snow_melt    ! Snow melt (kg H2O/m2/s)
%
% Output from Obukhov length calculation
%   fluxvar.ustar       ! Friction velocity (m/s)
%   fluxvar.tstar       ! Temperature scale (K)
%   fluxvar.qstar       ! Water vapor scale (mol/mol)
%   fluxvar.obu         ! Obukhov length (m)
%   fluxvar.z0m         ! RSL only: Roughness length for momentum (m)
%   fluxvar.z0c         ! RSL only: Roughness length for scalars (m)
%   fluxvar.disp        ! RSL only: Displacement height (m)
% ------------------------------------------------------

% --- Calculate the Obukhov length

% Calculate the Obukhov length (obu) for the current surface temperature
% and surface vapor pressure using Monin-Obukhov similarity theory or
% Harman & Finnigan (2007, 2008) roughness sublayer (RSL) theory. Use the
% functions "most" or "rsl" to iterate obu until the change in obu is less
% than tol.

obu0 = 100;                        % Initial estimate for Obukhov length (m)
obu1 = -100;                       % Initial estimate for Obukhov length (m)
tol = 0.01;                        % Accuracy tolerance for Obukhov length (m)

switch fluxvar.profiles
   case 'MOST'                     % Use Monin-Obukhov similarity theory
   func_name = 'most';             % The function name is "most", in the file most.m
   case 'RSL'                      % Use canopy coupling with roughness sublayer theory
   func_name = 'rsl';              % The function name is "rsl", in the file rsl.m
end

% Solve for the Obukhov length

[fluxvar, oburoot] = hybrid_root (func_name, physcon, forcvar, surfvar, fluxvar, obu0, obu1, tol);

% Uncomment this line to use MOST or RSL for neutral conditions
% [fluxvar, dummy] = most (physcon, forcvar, surfvar, fluxvar, -inf);
% [fluxvar, dummy] = rsl (physcon, forcvar, surfvar, fluxvar, -inf);

% --- Aerodynamic conductances for momentum (gam) and scalars (gac) (mol/m2/s)

fluxvar.gam = forcvar.rhomol * fluxvar.ustar * fluxvar.ustar / forcvar.uref;
fluxvar.gac = forcvar.rhomol * fluxvar.ustar * fluxvar.tstar / (forcvar.thref - fluxvar.tsrf);

% --- Calculate soil temperatures

% Save current temperatures

for i = 1:soilvar.nsoi
   tsoi0(i) = soilvar.tsoi(i);
end

% Saturation vapor pressure (Pa) and temperature derivative (Pa/K)

[esat, desat] = satvap (fluxvar.tsrf-physcon.tfrz);

% Latent heat of vaporization or sublimation (J/mol)

if (bucket.snow_water > 0)
   lambda_val = physcon.hsub;              % Latent heat of sublimation (J/kg)
else
   lambda_val = physcon.hvap;              % Latent heat of vaporization (J/kg)
end
lambda_val = lambda_val * physcon.mmh2o;   % J/kg -> J/mol

% Surface conductance for water vapor (mol/m2/s)

gw = 1 / (1/surfvar.gcan + 1/fluxvar.gac);

% Soil wetness factor for evapotranspiration

switch fluxvar.bucket
   case 'no_bucket'                % No bucket hydrology
   bucket.soil_beta = 1;
   case 'use_bucket'               % Use bucket hydrology
   bucket.soil_beta = min(bucket.soil_water/(bucket.soil_beta_max*bucket.soil_water_max), 1);
end

% Emitted longwave radiation (W/m2) and temperature derivative (W/m2/K)

fluxvar.lwrad = surfvar.emiss * physcon.sigma * fluxvar.tsrf^4;
dlwrad = 4 * surfvar.emiss * physcon.sigma * fluxvar.tsrf^3;

% Sensible heat flux (W/m2) and temperature derivative (W/m2/K)

fluxvar.shflx = forcvar.cpair * (fluxvar.tsrf - forcvar.thref) * fluxvar.gac;
dshflx = forcvar.cpair * fluxvar.gac;

% Latent heat flux (W/m2) and temperature derivative (W/m2/K)

fluxvar.lhflx = lambda_val / forcvar.pref * (esat - forcvar.eref) * gw * bucket.soil_beta;
dlhflx = lambda_val / forcvar.pref * desat * gw * bucket.soil_beta;

% Net energy flux into soil (W/m2) and temperature derivative (W/m2/K)

f0 = fluxvar.qa - fluxvar.lwrad - fluxvar.shflx - fluxvar.lhflx;
df0 = -dlwrad - dshflx - dlhflx;

% Update soil temperatures

[soilvar, fluxvar, bucket] = soil_temperature (physcon, soilvar, fluxvar, bucket, dt, f0, df0);

% --- Update surface fluxes for the change in surface temperature

dtsrf = soilvar.tsoi(1) - tsoi0(1);
fluxvar.lwrad = fluxvar.lwrad + dlwrad * dtsrf;
fluxvar.shflx = fluxvar.shflx + dshflx * dtsrf;
fluxvar.lhflx = fluxvar.lhflx + dlhflx * dtsrf;
fluxvar.gsoi = f0 + df0 * dtsrf;

% --- Adjust heat flux into soil for snow melt

fluxvar.gsoi = fluxvar.gsoi - fluxvar.gsno;

% --- Net radiation

fluxvar.rnet = fluxvar.qa - fluxvar.lwrad;

% --- Error check

err = fluxvar.rnet - fluxvar.shflx - fluxvar.lhflx - fluxvar.gsoi - fluxvar.gsno;
if (abs(err) > 1e-06)
   fprintf('err = %15.3f\n',err)
   fprintf('qa = %15.3f\n',fluxvar.qa)
   fprintf('lwrad = %15.3f\n',fluxvar.lwrad)
   fprintf('sh = %15.3f\n',fluxvar.shflx)
   fprintf('lh = %15.3f\n',fluxvar.lhflx)
   fprintf('gsoi = %15.3f\n',fluxvar.gsoi)
   fprintf('gsno = %15.3f\n',fluxvar.gsno)
   error ('surface temperature error')
end

% --- Evapotranspiration (mol H2O/m2/s)

fluxvar.etflx = fluxvar.lhflx / lambda_val;

% --- Surface vapor pressure is diagnosed from evaporative flux

fluxvar.esrf = (forcvar.eref / forcvar.pref) + fluxvar.etflx / fluxvar.gac; % mol/mol
fluxvar.esrf = fluxvar.esrf * forcvar.pref;                                 % mol/mol -> Pa

% --- Phase change for soil layers undergoing freezing of thawing

switch soilvar.method

   case 'apparent-heat-capacity'

   % No explicit phase change energy flux. This is included in the heat capacity.

   soilvar.hfsoi = 0;

   case 'excess-heat'

   % Adjust temperatures for phase change. Freeze or melt ice using energy
   % excess or deficit needed to change temperature to the freezing point.
   % The variable hfsoi is returned as the energy flux from phase change (W/m2).

   [soilvar] = phase_change (physcon, soilvar, dt);

end

% --- Check for energy conservation

% Sum change in energy (W/m2)

edif = 0;
for i = 1:soilvar.nsoi
   edif = edif + soilvar.cv(i) * soilvar.dz(i) * (soilvar.tsoi(i) - tsoi0(i)) / dt;
end

% Error check

err = edif - fluxvar.gsoi - soilvar.hfsoi;
if (abs(err) > 1e-03)
   error ('Soil temperature energy conservation error')
end

% --- Surface temperature is the first soil layer

fluxvar.tsrf = soilvar.tsoi(1);
