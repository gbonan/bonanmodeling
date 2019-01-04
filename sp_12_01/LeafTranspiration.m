function [flux] = LeafTranspiration (physcon, atmos, flux)

% Leaf transpiration flux

% ------------------------------------------------------
% Input
%   physcon.tfrz     ! Freezing point of water (K)
%   physcon.mmh2o    ! Molecular mass of water (kg/mol)
%   atmos.patm       ! Atmospheric pressure (Pa)
%   atmos.tair       ! Air temperature (K)
%   atmos.eair       ! Vapor pressure of air (Pa)
%   flux.gbv         ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
%   flux.gs          ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
%   flux.tleaf       ! Leaf temperature (K)
%
% Output
%   flux.lhflx       ! Leaf latent heat flux (W/m2 leaf)
%   flux.etflx       ! Leaf transpiration flux (mol H2O/m2 leaf/s)
% ------------------------------------------------------

% Latent heat of vaporization (J/mol)

[lambda_val] = latvap ((atmos.tair-physcon.tfrz), physcon.mmh2o);

% Saturation vapor pressure (Pa) and temperature derivative (Pa/K)

[esat, desat] = satvap ((flux.tleaf-physcon.tfrz));

% Leaf conductance for water vapor (mol H2O/m2/s)

gleaf = flux.gs * flux.gbv / (flux.gs + flux.gbv);

% Latent heat flux (W/m2)

flux.lhflx = lambda_val / atmos.patm * (esat - atmos.eair) * gleaf;

% Water vapor flux: W/m2 -> mol H2O/m2/s

flux.etflx = flux.lhflx / lambda_val;
