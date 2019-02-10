function [fluxvar] = canopy_turbulence (dt, physcon, forcvar, surfvar, leafvar, soilvar, fluxvar)

% Canopy turbulence, aeorodynamic conductances, and wind/temperature/water vapor
% profiles using above- and within-canopy coupling with a roughness sublayer
% (RSL) parameterization 

% ------------------------------------------------------------------------------------
% Input
%   dt                  ! Model time step (s)
%   physcon.vkc         ! von Karman constant
%   physcon.grav        ! Gravitational acceleration (m/s2)
%   physcon.mmh2o       ! Molecular mass of water (kg/mol)
%   physcon.tfrz        ! Freezing point of water (K)
%   physcon.hvap        ! Latent heat of evaporation (J/kg)
%
%   forcvar.zref        ! Reference height (m)
%   forcvar.uref        ! Wind speed at reference height (m/s)
%   forcvar.thref       ! Potential temperature at reference height (K)
%   forcvar.thvref      ! Virtual potential temperature at reference height (K)
%   forcvar.qref        ! Water vapor at reference height (mol/mol)
%   forcvar.pref        ! Atmospheric pressure (Pa)
%   forcvar.mmair       ! Molecular mass of air at reference height (kg/mol)
%   forcvar.cpair       ! Specific heat of air at constant pressure (J/mol/K)
%   forcvar.rhomol      ! Molar density (mol/m3)
%
%   surfvar.p           ! Index for grid point to process
%   surfvar.hc          ! Canopy height (m)
%   surfvar.pai         ! Canopy plant area index (m2/m2)
%   surfvar.nlev        ! Index for top level
%   surfvar.ntop        ! Index for top leaf layer
%   surfvar.nsoi        ! First canopy layer is soil
%   surfvar.zs          ! Canopy height for scalar concentration and source (m)
%   surfvar.zw          ! Canopy height at layer interfaces (m)
%   surfvar.dpai        ! Layer plant area index (m2/m2)
%   surfvar.fwet        ! Fraction of plant area index that is wet
%   surfvar.fdry        ! Fraction of plant area index that is green and dry
%   surfvar.fracsun     ! Sunlit fraction of canopy layer
%   surfvar.fracsha     ! Shaded fraction of canopy layer
%
%   leafvar.nleaf       ! Number of leaf types (sunlit and shaded)
%   leafvar.isun        ! Sunlit leaf index
%   leafvar.isha        ! Shaded leaf index
%   leafvar.gbh         ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
%   leafvar.gbv         ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
%   leafvar.gs          ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
%   leafvar.cpleaf      ! Leaf heat capacity (J/m2 leaf/K)
%   leafvar.rnleaf      ! Leaf net radiation (W/m2 leaf)
%
%   soilvar.tk          ! Soil thermal conductivity (W/m/K)
%   soilvar.dz          ! Soil layer depth (m)
%   soilvar.tsoi        ! Soil temperature (K)
%   soilvar.resis       ! Soil evaporative resistance (s/m)
%   soilvar.rhg         ! Relative humidity of airspace at soil surface (fraction)
%
%   fluxvar.rnsoi       ! Net radiation at ground (W/m2)
%   fluxvar.tair_old    ! Air temperature profile for previous timestep (K)
%   fluxvar.qair_old    ! Water vapor profile for previous timestep (mol/mol)
%   fluxvar.tveg_old    ! Vegetation temperature profile for previous timestep (K)
%
% Input/output
%   fluxvar.taf         ! Air temperature at canopy top (K)
%   fluxvar.qaf         ! Water vapor at canopy top (mol/mol)
%
% Output (specific to canopy_turbulence)
%   fluxvar.Lc          ! Canopy density length scale (m)
%   fluxvar.wind        ! Wind speed profile (m/s)
%   fluxvar.uaf         ! Wind speed at canopy top (m/s)
%   fluxvar.ga_prof     ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
%
% Output (from obukhov_function)
%   fluxvar.c1m         ! Roughness sublayer c1 parameter for momentum (dimensionless)
%   fluxvar.c1c         ! Roughness sublayer c1 parameter for scalars (dimensionless)
%   fluxvar.c2          ! Roughness sublayer depth scale multiplier (dimensionless)
%   fluxvar.disp        ! Displacement height (m)
%   fluxvar.beta        ! u* / u(hc)
%   fluxvar.PrSc        ! Prandtl (Schmidt) number at canopy top
%   fluxvar.ustar       ! Friction velocity (m/s)
%   fluxvar.tstar       ! Temperature scale (K)
%   fluxvar.qstar       ! Water vapor scale (mol/mol)
%   fluxvar.gac         ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
%   fluxvar.obu_ustar   ! Obukhov length used for u* (m)
%   fluxvar.obu         ! Value for Obukhov length (m)
%
% Output (from scalar_profile)
%   fluxvar.tair        ! Air temperature profile (K)
%   fluxvar.qair        ! Water vapor profile (mol/mol)
%   fluxvar.tveg        ! Vegetation temperature profile (K)
%   fluxvar.shsoi       ! Ground sensible heat flux, ground (W/m2)
%   fluxvar.etsoi       ! Ground evaporation flux (mol H2O/m2/s)
%   fluxvar.gsoi        ! Soil heat flux (W/m2)
%   fluxvar.tg          ! Soil surface temperature (K)
%   fluxvar.shveg       ! Leaf sensible heat flux (W/m2 leaf)
%   fluxvar.etveg       ! Leaf evapotranspiration flux (mol H2O/m2 leaf/s)
%   fluxvar.stveg       ! Leaf storage heat flux (W/m2 leaf)
%   fluxvar.shair       ! Canopy air sensible heat flux (W/m2)
%   fluxvar.etair       ! Canopy air water vapor flux (mol H2O/m2/s)
%   fluxvar.stair       ! Canopy air storage heat flux (W/m2)
% ------------------------------------------------------------------------------------

% --- Index for grid point to process

p = surfvar.p;

% --- Leaf drag coefficient (dimensionless)

cd = 0.25;

% --- Canopy density length scale (m)

fluxvar.Lc(p) = surfvar.hc(p) / (cd * surfvar.pai(p));

% --- Calculate the Obukhov length

% Calculate the Obukhov length (obu) for the current surface temperature
% and surface vapor pressure using Harman & Finnigan (2007, 2008) roughness 
% sublayer (RSL) theory. Use the function "obukhov_function" to iterate obu
% until the change in obu is less than tol.

obu0 = 100;                     % Initial estimate for Obukhov length (m)
obu1 = -100;                    % Initial estimate for Obukhov length (m)
tol = 0.01;                     % Accuracy tolerance for Obukhov length (m)
func_name = 'obukhov_function'; % The function is the file obukhov_function.m

% Solve for the Obukhov length. Do not use final returned value for obu.
% Instead, use the value used to calculate u*

[fluxvar, oburoot] = hybrid_root (func_name, physcon, forcvar, surfvar, fluxvar, obu0, obu1, tol);
fluxvar.obu(p) = fluxvar.obu_ustar(p);

% --- Wind profile (m/s)

% Above-canopy wind speed: defined at zs

h_minus_d = surfvar.hc(p) - fluxvar.disp(p);
[psi_m_hc] = psi_m_monin_obukhov (h_minus_d / fluxvar.obu(p));
[psi_m_rsl_hc] = psi_m_rsl (h_minus_d, h_minus_d, fluxvar.obu(p), fluxvar.c1m(p), fluxvar.c2);

for ic = surfvar.ntop(p)+1:surfvar.nlev(p)
   z_minus_d = surfvar.zs(p,ic) - fluxvar.disp(p);
   [psi_m_zs] = psi_m_monin_obukhov (z_minus_d / fluxvar.obu(p));
   [psi_m_rsl_zs] = psi_m_rsl (z_minus_d, h_minus_d, fluxvar.obu(p), fluxvar.c1m(p), fluxvar.c2);
   psim = -psi_m_zs + psi_m_hc + psi_m_rsl_zs - psi_m_rsl_hc + physcon.vkc / fluxvar.beta(p);
   fluxvar.wind(p,ic) = fluxvar.ustar(p) / physcon.vkc * (log(z_minus_d/h_minus_d) + psim);
end

% Wind speed at top of canopy

fluxvar.uaf(p) = fluxvar.ustar(p) / fluxvar.beta(p);

% Within-canopy wind speed: defined at zs. Limit to > 0.1 m/s

lm = 2 * fluxvar.beta(p)^3 * fluxvar.Lc(p);
lm_over_beta = lm / fluxvar.beta(p);
for ic = surfvar.nsoi(p)+1:surfvar.ntop(p)
   fluxvar.wind(p,ic) = fluxvar.uaf(p) * exp((surfvar.zs(p,ic) - surfvar.hc(p)) / lm_over_beta);
   fluxvar.wind(p,ic) = max(fluxvar.wind(p,ic), 0.1);
end

% Wind speed at ground

fluxvar.wind(p,surfvar.nsoi(p)) = 0;

% --- Aerodynamic conductances (mol/m2/s) between zs(i) and zs(i+1)

% Above-canopy conductances

h_minus_d = surfvar.hc(p) - fluxvar.disp(p);
[psi_c_hc] = psi_c_monin_obukhov (h_minus_d / fluxvar.obu(p));
[psi_c_rsl_hc] = psi_c_rsl (h_minus_d, h_minus_d, fluxvar.obu(p), fluxvar.c1c(p), fluxvar.c2);

for ic = surfvar.ntop(p)+1:surfvar.nlev(p)-1

   % Lower height zs(i)

   z_minus_d = surfvar.zs(p,ic) - fluxvar.disp(p);
   [psi_c_zs] = psi_c_monin_obukhov (z_minus_d / fluxvar.obu(p));
   [psi_c_rsl_zs] = psi_c_rsl (z_minus_d, h_minus_d, fluxvar.obu(p), fluxvar.c1c(p), fluxvar.c2);
   psic1 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc;

   % Upper height zs(i+1)

   z_minus_d = surfvar.zs(p,ic+1) - fluxvar.disp(p);
   [psi_c_zs] = psi_c_monin_obukhov (z_minus_d / fluxvar.obu(p));
   [psi_c_rsl_zs] = psi_c_rsl (z_minus_d, h_minus_d, fluxvar.obu(p), fluxvar.c1c(p), fluxvar.c2);
   psic2 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc;

   % Conductance. Note that psic = psic2 - psic1 is equivalent to:
   % -psi_c_z2 + psi_c_z1 + psi_c_rsl_z2 - psi_c_rsl_z1

   psic = psic2 - psic1;
   zlog = log((surfvar.zs(p,ic+1)-fluxvar.disp(p)) / (surfvar.zs(p,ic)-fluxvar.disp(p)));
   fluxvar.ga_prof(p,ic) = forcvar.rhomol(p) * physcon.vkc * fluxvar.ustar(p) / (zlog + psic);
end

% Special case for the top layer to the reference height

ic = surfvar.nlev(p);
z_minus_d = surfvar.zs(p,ic) - fluxvar.disp(p);
[psi_c_zs] = psi_c_monin_obukhov (z_minus_d / fluxvar.obu(p));
[psi_c_rsl_zs] = psi_c_rsl (z_minus_d, h_minus_d, fluxvar.obu(p), fluxvar.c1c(p), fluxvar.c2);
psic1 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc;

z_minus_d = forcvar.zref(p) - fluxvar.disp(p);
[psi_c_zs] = psi_c_monin_obukhov (z_minus_d / fluxvar.obu(p));
[psi_c_rsl_zs] = psi_c_rsl (z_minus_d, h_minus_d, fluxvar.obu(p), fluxvar.c1c(p), fluxvar.c2);
psic2 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc;

psic = psic2 - psic1;
zlog = log((forcvar.zref(p)-fluxvar.disp(p)) / (surfvar.zs(p,ic)-fluxvar.disp(p)));
fluxvar.ga_prof(p,ic) = forcvar.rhomol(p) * physcon.vkc * fluxvar.ustar(p) / (zlog + psic);

% Within-canopy aerodynamic conductances

for ic = surfvar.nsoi(p)+1:surfvar.ntop(p)-1
   zl = surfvar.zs(p,ic) - surfvar.hc(p);
   zu = surfvar.zs(p,ic+1) - surfvar.hc(p);
   res = fluxvar.PrSc(p) / (fluxvar.beta(p) * fluxvar.ustar(p)) * (exp(-zl/lm_over_beta) - exp(-zu/lm_over_beta));
   fluxvar.ga_prof(p,ic) = forcvar.rhomol(p) / res;
end

% Special case for top canopy layer: conductance from zs to hc

ic = surfvar.ntop(p);
zl = surfvar.zs(p,ic) - surfvar.hc(p);
zu = surfvar.hc(p) - surfvar.hc(p);
res = fluxvar.PrSc(p) / (fluxvar.beta(p) * fluxvar.ustar(p)) * (exp(-zl/lm_over_beta) - exp(-zu/lm_over_beta));
ga_below_hc = forcvar.rhomol(p) / res;

% Now include additional conductance from hc to first atmospheric layer

z_minus_d = surfvar.hc(p) - fluxvar.disp(p);
[psi_c_zs] = psi_c_monin_obukhov (z_minus_d / fluxvar.obu(p));
[psi_c_rsl_zs] = psi_c_rsl (z_minus_d, h_minus_d, fluxvar.obu(p), fluxvar.c1c(p), fluxvar.c2);
psic1 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc;

z_minus_d = surfvar.zs(p,ic+1) - fluxvar.disp(p);
[psi_c_zs] = psi_c_monin_obukhov (z_minus_d / fluxvar.obu(p));
[psi_c_rsl_zs] = psi_c_rsl (z_minus_d, h_minus_d, fluxvar.obu(p), fluxvar.c1c(p), fluxvar.c2);
psic2 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc;

psic = psic2 - psic1;
zlog = log((surfvar.zs(p,ic+1)-fluxvar.disp(p)) / (surfvar.hc(p)-fluxvar.disp(p)));
ga_above_hc = forcvar.rhomol(p) * physcon.vkc * fluxvar.ustar(p) / (zlog + psic);
fluxvar.ga_prof(p,ic) = 1 / (1 / ga_below_hc + 1 / ga_above_hc);

% Make sure above-canopy aerodynamic resistances sum to 1/gac

sumres = 1 / ga_above_hc;
for ic = surfvar.ntop(p)+1:surfvar.nlev(p)
   sumres = sumres + 1 / fluxvar.ga_prof(p,ic);
end

if (abs(1/sumres - fluxvar.gac(p)) > 1e-06)
   error('canopy_turbulence: above-canopy aerodynamic conductance error')
end

% Aerodynamic conductance at ground

ic = surfvar.nsoi(p);
ic_plus_one = surfvar.nsoi(p)+1;

z0m_g = 0.01;                         % Roughness length of ground (m)
z0c_g = 0.1 * z0m_g;                  % Roughness length for scalars
zlog_m = log(surfvar.zs(p,ic_plus_one)/z0m_g);
zlog_c = log(surfvar.zs(p,ic_plus_one)/z0c_g);
ustar_g = fluxvar.wind(p,ic_plus_one) * physcon.vkc / zlog_m;
ustar_g = max(ustar_g, 0.01);
res = zlog_c / (physcon.vkc * ustar_g);
fluxvar.ga_prof(p,ic) = forcvar.rhomol(p) / res;

% --- Limit resistances to < 500 s/m

for ic = surfvar.nsoi(p):surfvar.nlev(p)
   res = min (forcvar.rhomol(p)/fluxvar.ga_prof(p,ic), 500);
   fluxvar.ga_prof(p,ic) = forcvar.rhomol(p) / res;
end

% --- Calculate within-canopy scalar profiles for temperature and vapor pressure

[fluxvar] = scalar_profile (dt, p, physcon, forcvar, surfvar, leafvar, soilvar, fluxvar);

% --- Temperature and water vapor at top of canopy

fluxvar.taf(p) = fluxvar.tair(p,surfvar.ntop(p));
fluxvar.qaf(p) = fluxvar.qair(p,surfvar.ntop(p));
