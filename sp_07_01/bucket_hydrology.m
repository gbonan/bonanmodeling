function [bucket] = bucket_hydrology (physcon, forcvar, fluxvar, bucket, dt)

% Bucket model hydrology

% ------------------------------------------------------
% Input
%   dt                    ! Time step (s)
%   physcon.mmh2o         ! Molecular mass of water (kg/mol)
%   forcvar.rain          ! Rainfall (kg H2O/m2/s)
%   forcvar.snow          ! Snowfall (kg H2O/m2/s)
%   fluxvar.etflx         ! Evapotranspiration (mol H2O/m2/s)
%   bucket.snow_melt      ! Snow melt (kg H2O/m2/s)
%   bucket.soil_water_max ! Maximum soil water (kg H2O/m2)
%
% Input/output
%   bucket.snow_water     ! Snow water (kg H2O/m2)
%   bucket.soil_water     ! Soil water (kg H2O/m2)
%
% Output
%   bucket.runoff         ! Runoff (kg H2O/m2/s)
% ------------------------------------------------------

% --- Save current water for conservation check

snow0 = bucket.snow_water;
soil0 = bucket.soil_water;
tot0 = snow0 + soil0;

% --- Update snow and soil water for snow melt

bucket.snow_water = bucket.snow_water - bucket.snow_melt * dt;
bucket.soil_water = bucket.soil_water + bucket.snow_melt * dt;

% --- Update snow and soil water for precipitation

bucket.soil_water = bucket.soil_water + forcvar.rain * dt;
bucket.snow_water = bucket.snow_water + forcvar.snow * dt;

% --- Update snow and soil water for evaporative flux

% Evaporative loss (kg H2O/m2/s)

evap = fluxvar.etflx * physcon.mmh2o;

% Apply evaporative flux to snow, but if more snow sublimates than
% is present take the extra water from the soil as evaporation

subl = min (evap, bucket.snow_water/dt);
evap = evap - subl;

% Update snow and soil

bucket.snow_water = bucket.snow_water - subl * dt;
bucket.soil_water = bucket.soil_water - evap * dt;

% --- Runoff is soil water in excess of a maximum capacity

bucket.runoff = max ((bucket.soil_water - bucket.soil_water_max)/dt, 0);
bucket.soil_water = bucket.soil_water - bucket.runoff * dt;

% --- Conservation check

% Snow water

delta = bucket.snow_water - snow0;
err = (forcvar.snow - bucket.snow_melt - subl) * dt - delta;
if (abs(err) > 1e-03)
   error ('Bucket model snow conservation error')
end

% Soil water

delta = bucket.soil_water - soil0;
err = (forcvar.rain + bucket.snow_melt - evap - bucket.runoff) * dt - delta;
if (abs(err) > 1e-03)
   error ('Bucket model soil conservation error')
end

% Total water

delta = (bucket.snow_water + bucket.soil_water) - tot0;
err = (forcvar.rain + forcvar.snow - fluxvar.etflx * physcon.mmh2o - bucket.runoff) * dt - delta;
if (abs(err) > 1e-03)
   error ('Bucket model total conservation error')
end
