function [soilvar, fluxvar, bucket] = soil_temperature (physcon, soilvar, fluxvar, bucket, dt, f0, df0)

% Use an implicit formulation with the surface boundary condition specified
% as the soil heat flux (W/m2; positive into the soil) to solve for soil
% temperatures at time n+1.
%
% Calculate soil temperatures as:
%
%      dT   d     dT 
%   cv -- = -- (k --)
%      dt   dz    dz 
%
% where: T = temperature (K)
%        t = time (s)
%        z = depth (m)
%        cv = volumetric heat capacity (J/m3/K)
%        k = thermal conductivity (W/m/K)
%
% Set up a tridiagonal system of equations to solve for dT at time n+1, 
% where the temperature equation for layer i is:
%
%   d_i = a_i*[dT_i-1] n+1 + b_i*[dT_i] n+1 + c_i*[dT_i+1] n+1
%
% The energy flux into the soil is specified as a linear function of
% the temperature of the first soil layer T(1) using the equation:
%
%   gsoi = f0([T_1] n) + df0 / dT * ([T_1] n+1 - [T_1] n)
%
% where f0 is the surface energy balance evaluated with T(1) at time n
% and df0/dT is the temperature derivative of f0.
%
% If snow is present and the surface temperature is above freezing, reset
% the temperature to freezing and use the residual energy to melt snow.

% ------------------------------------------------------
% Input
%   dt                      ! Time step (s)
%   f0                      ! Energy flux into soil (W/m2)
%   df0                     ! Temperature derivative of f0 (W/m2/K)
%   physcon.hfus            ! Heat of fusion for water at 0 C (J/kg)
%   physcon.tfrz            ! Freezing point of water (K)
%   soilvar.nsoi            ! Number of soil layers
%   soilvar.z               ! Soil depth (m)
%   soilvar.z_plus_onehalf  ! Soil depth (m) at i+1/2 interface between layers i and i+1
%   soilvar.dz              ! Soil layer thickness (m)
%   soilvar.dz_plus_onehalf ! Thickness (m) between between i and i+1
%   soilvar.tk              ! Thermal conductivity (W/m/K)
%   soilvar.cv              ! Heat capacity (J/m3/K)
%   bucket.snow_water       ! Snow water (kg H2O/m2)
%
% Input/output
%   soilvar.tsoi            ! Soil temperature (K)
%
% Output
%   fluxvar.gsno            ! Snow melt energy flux (W/m2)
%   bucket.snow_melt        ! Snow melt (kg H2O/m2/s)
% ------------------------------------------------------

% --- Save current soil temperature

for i = 1:soilvar.nsoi
   tsoi0(i) = soilvar.tsoi(i);
end

% --- Thermal conductivity at interface (W/m/K)

for i = 1:soilvar.nsoi-1
   tk_plus_onehalf(i) = soilvar.tk(i) * soilvar.tk(i+1) * (soilvar.z(i)-soilvar.z(i+1)) / ...
   (soilvar.tk(i)*(soilvar.z_plus_onehalf(i)-soilvar.z(i+1)) + soilvar.tk(i+1)*(soilvar.z(i)-soilvar.z_plus_onehalf(i)));
end

% --- Set up tridiagonal matrix

% Top soil layer

i = 1;
a(i) = 0;
c(i) = -tk_plus_onehalf(i) / soilvar.dz_plus_onehalf(i);
b(i) = soilvar.cv(i) * soilvar.dz(i) / dt - c(i) - df0;
d(i) = -tk_plus_onehalf(i) * (soilvar.tsoi(i) - soilvar.tsoi(i+1)) / soilvar.dz_plus_onehalf(i) + f0;

% Layers 2 to nsoi-1

for i = 2:soilvar.nsoi-1
   a(i) = -tk_plus_onehalf(i-1) / soilvar.dz_plus_onehalf(i-1);
   c(i) = -tk_plus_onehalf(i) / soilvar.dz_plus_onehalf(i);
   b(i) = soilvar.cv(i) * soilvar.dz(i) / dt - a(i) - c(i);
   d(i) = tk_plus_onehalf(i-1) * (soilvar.tsoi(i-1) - soilvar.tsoi(i)) / soilvar.dz_plus_onehalf(i-1) ...
        - tk_plus_onehalf(i) * (soilvar.tsoi(i) - soilvar.tsoi(i+1)) / soilvar.dz_plus_onehalf(i);
end

% Bottom soil layer

i = soilvar.nsoi;
a(i) = -tk_plus_onehalf(i-1) / soilvar.dz_plus_onehalf(i-1);
c(i) = 0;
b(i) = soilvar.cv(i) * soilvar.dz(i) / dt - a(i);
d(i) = tk_plus_onehalf(i-1) * (soilvar.tsoi(i-1) - soilvar.tsoi(i)) / soilvar.dz_plus_onehalf(i-1);

% --- Begin tridiagonal solution: forward sweep for layers N to 1

% Bottom soil layer

i = soilvar.nsoi;
e(i) = a(i) / b(i);
f(i) = d(i) / b(i);

% Layers nsoi-1 to 2

for i = soilvar.nsoi-1: -1: 2
   den = b(i) - c(i) * e(i+1);
   e(i) = a(i) / den;
   f(i) = (d(i) - c(i) * f(i+1)) / den;
end

% Complete the tridiagonal solution to get the temperature of the top soil layer

i = 1;
num = d(i) - c(i) * f(i+1);
den = b(i) - c(i) * e(i+1);
tsoi_test = soilvar.tsoi(i) + num / den;

% --- Surface temperature with adjustment for snow melt

% Potential melt rate based on temperature above freezing

potential_snow_melt = max(0, (tsoi_test - physcon.tfrz) * den / physcon.hfus);

% Maximum melt rate is the amount of snow that is present

maximum_snow_melt = bucket.snow_water / dt;

% Cannot melt more snow than is present

bucket.snow_melt = min(maximum_snow_melt, potential_snow_melt);

% Energy flux for snow melt

fluxvar.gsno = bucket.snow_melt * physcon.hfus;

% Update temperature - If there is no snow melt, tsoi(1) = tsoi_test (as above).
% While snow is melting at the potential rate, tsoi(1) = tfrz. If snow melt is
% less than the potential rate, tsoi(1) > tfrz and < tsoi_test.

soilvar.tsoi(i) = soilvar.tsoi(i) + (num - fluxvar.gsno) / den;
dtsoi(i) = soilvar.tsoi(i) - tsoi0(i);

% --- Now complete the tridiagonal solution for layers 2 to N

for i = 2:soilvar.nsoi
   dtsoi(i) = f(i) - e(i) * dtsoi(i-1);
   soilvar.tsoi(i) = soilvar.tsoi(i) + dtsoi(i);
end
