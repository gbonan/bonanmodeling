function [flux] = LeafTemperature (physcon, atmos, leaf, flux)

% Leaf temperature and energy fluxes

% ------------------------------------------------------
% Input
%   physcon.tfrz     ! Freezing point of water (K)
%   physcon.mmh2o    ! Molecular mass of water (kg/mol)
%   physcon.sigma    ! Stefan-Boltzmann constant (W/m2/K4)
%   atmos.patm       ! Atmospheric pressure (Pa)
%   atmos.cpair      ! Specific heat of air at constant pressure (J/mol/K)
%   atmos.tair       ! Air temperature (K)
%   atmos.eair       ! Vapor pressure of air (Pa)
%   leaf.emiss       ! Leaf emissivity
%   flux.gbh         ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
%   flux.gbv         ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
%   flux.gs          ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
%   flux.qa          ! Leaf radiative forcing (W/m2 leaf)
%
% Input/ouput
%   flux.tleaf       ! Leaf temperature (K)
%
% Output
%   flux.rnet        ! Leaf net radiation (W/m2 leaf)
%   flux.lwrad       ! Longwave radiation emitted from leaf (W/m2 leaf)
%   flux.shflx       ! Leaf sensible heat flux (W/m2 leaf)
%   flux.lhflx       ! Leaf latent heat flux (W/m2 leaf)
%   flux.etflx       ! Leaf transpiration flux (mol H2O/m2 leaf/s)
% ------------------------------------------------------

% --- Latent heat of vaporization (J/mol)

[lambda_val] = latvap ((atmos.tair-physcon.tfrz), physcon.mmh2o);

% --- Newton-Raphson iteration until leaf energy balance is less than
% f0_max or to niter_max iterations

niter = 0;         % Number of iterations
f0 = 1e36;         % Leaf energy balance (W/m2)

niter_max = 100;   % Maximum number of iterations
f0_max = 1e-06;    % Maximum energy imbalance (W/m2)

while (niter <= niter_max & abs(f0) > f0_max)

   % Increment iteration counter

   niter = niter + 1;

   % Saturation vapor pressure (Pa) and temperature derivative (Pa/K)

   [esat, desat] = satvap ((flux.tleaf-physcon.tfrz));

   % Leaf conductance for water vapor (mol H2O/m2/s)

   gleaf = flux.gs * flux.gbv / (flux.gs + flux.gbv);

   % Emitted longwave radiation (W/m2) and temperature derivative (W/m2/K)

   flux.lwrad = 2 * leaf.emiss * physcon.sigma * flux.tleaf^4;
   dlwrad = 8 * leaf.emiss * physcon.sigma * flux.tleaf^3;

   % Sensible heat flux (W/m2) and temperature derivative (W/m2/K)

   flux.shflx = 2 * atmos.cpair * (flux.tleaf - atmos.tair) * flux.gbh;
   dshflx = 2 * atmos.cpair * flux.gbh;

   % Latent heat flux (W/m2) and temperature derivative (W/m2/K)

   flux.lhflx = lambda_val / atmos.patm * (esat - atmos.eair) * gleaf;
   dlhflx = lambda_val / atmos.patm * desat * gleaf;

   % Energy balance (W/m2) and temperature derivative (W/m2/K)

   f0 = flux.qa - flux.lwrad - flux.shflx - flux.lhflx;
   df0 = -dlwrad - dshflx - dlhflx;

   % Change in leaf temperature

   dtleaf = -f0 / df0;

   % Update leaf temperature

   flux.tleaf = flux.tleaf + dtleaf;

end

% --- Net radiation

flux.rnet = flux.qa - flux.lwrad;

% --- Error check

err = flux.rnet - flux.shflx - flux.lhflx;
if (abs(err) > f0_max)
   fprintf('err = %15.3f\n',err)
   fprintf('qa = %15.3f\n',flux.qa)
   fprintf('lwrad = %15.3f\n',flux.lwrad)
   fprintf('sh = %15.3f\n',flux.shflx)
   fprintf('lh = %15.3f\n',flux.lhflx)
   error ('LeafTemperature error')
end

% Water vapor flux: W/m2 -> mol H2O/m2/s

flux.etflx = flux.lhflx / lambda_val;
