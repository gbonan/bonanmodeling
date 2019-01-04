function [flux] = SoilResistance (physcon, leaf, rootvar, soil, flux)

% Calculate soil hydraulic resistance, weighted soil water potential and
% water uptake from each soil layer

% ----------------------------------------------------------------------------------
% Input
%   physcon.grav      ! Gravitational acceleration (m/s2)
%   physcon.denh2o    ! Density of liquid water (kg/m3)
%   physcon.mmh2o     ! Molecular mass of water (kg/mol)
%   leaf.minlwp       ! Minimum leaf water potential (MPa)
%   rootvar.biomass   ! Fine root biomass (g biomass / m2)
%   rootvar.radius    ! Fine root radius (m)
%   rootvar.density   ! Fine root density (g biomass / m3 root)
%   rootvar.resist    ! Hydraulic resistivity of root tissue (MPa.s.g/mmol H2O)
%   soil.nlevsoi      ! Number of soil layers
%   soil.h2osoi_vol   ! Soil layer volumetric water content (m3/m3)
%   soil.psi          ! Soil layer matric potential (mm)
%   soil.watsat       ! Soil layer volumetric water content at saturation (porosity)
%   soil.hksat        ! Soil layer hydraulic conductivity at saturation (mm H2O/s)
%   soil.bsw          ! Soil layer Clapp and Hornberger "b" parameter
%   soil.rootfr       ! Fraction of roots in each soil layer (-)
%   soil.dz           ! Soil layer thickness (m)
%   flux.lai          ! Canopy leaf area index (m2/m2)
% Output
%   flux.rsoil        ! Soil hydraulic resistance (MPa.s.m2/mmol H2O)
%   flux.psi_soil     ! Weighted soil water potential (MPa)
%   flux.et_loss      ! Fraction of total transpiration from each soil layer (-)
% ----------------------------------------------------------------------------------

% --- Head of pressure  (MPa/m)

head = physcon.denh2o * physcon.grav * 1e-06;

% --- Root cross-sectional area (m2 root)

root_cross_sec_area = pi * rootvar.radius^2;

% --- Soil and root resistances for each layer

flux.rsoil = 0;
for j = 1:soil.nlevsoi

   % Hydraulic conductivity for each layer (mmol/m/s/MPa)

   s = max(min(soil.h2osoi_vol(j)/soil.watsat(j), 1), 0.01);
   hk = soil.hksat(j) * s^(2 * soil.bsw(j) + 3);     % mm/s
   hk = hk * 1e-03 / head;                           % mm/s -> m/s -> m2/s/MPa
   hk = hk * physcon.denh2o / physcon.mmh2o * 1000;  % m2/s/MPa -> mmol/m/s/MPa

   % Matric potential for each layer (MPa)

   psi_mpa(j) = soil.psi(j) * 1e-03 * head;          % mm -> m -> MPa

   % Root biomass density (g biomass / m3 soil)

   root_biomass_density = rootvar.biomass * soil.rootfr(j) / soil.dz(j);
   root_biomass_density = max(root_biomass_density, 1e-10);

   % Root length density (m root per m3 soil)

   root_length_density = root_biomass_density / (rootvar.density * root_cross_sec_area);

   % Distance between roots (m)

   root_dist = sqrt(1 / (root_length_density * pi));

   % Soil-to-root resistance (MPa.s.m2/mmol H2O)

   soilr1 = log(root_dist/rootvar.radius) / (2 * pi * root_length_density * soil.dz(j) * hk);

   % Root-to-stem resistance (MPa.s.m2/mmol H2O)

   soilr2 = rootvar.resist / (root_biomass_density * soil.dz(j));

   % Belowground resistance (MPa.s.m2/mmol H2O) 

   soilr = soilr1 + soilr2;

   % Total belowground resistance. First sum the conductances (1/soilr)
   % for each soil layer and then convert back to a resistance after the
   % summation

   flux.rsoil = flux.rsoil + 1 / soilr;

   % Maximum transpiration for each layer (mmol H2O/m2/s)

   evap(j) = (psi_mpa(j) - leaf.minlwp) / soilr;
   evap(j) = max (evap(j), 0);

end

% Belowground resistance: resistance = 1 / conductance

flux.rsoil = flux.lai / flux.rsoil;

% Weighted soil water potential (MPa)

flux.psi_soil = 0;
for j = 1:soil.nlevsoi
   flux.psi_soil = flux.psi_soil + psi_mpa(j) * evap(j);
end

totevap = sum(evap);      % Total maximum transpiration
if (totevap > 0)
   flux.psi_soil = flux.psi_soil / totevap;
else
   flux.psi_soil = leaf.minlwp;
end

% Fractional transpiration uptake from soil layers

for j = 1:soil.nlevsoi
   if (totevap > 0)
      flux.et_loss(j) = evap(j) / totevap;
   else
      flux.et_loss(j) = 1 / soil.nlevsoi;
   end
end
