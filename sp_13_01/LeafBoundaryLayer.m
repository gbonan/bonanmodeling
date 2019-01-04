function [flux] = LeafBoundaryLayer (physcon, atmos, leaf, flux)

% Leaf boundary layer conductances

% -------------------------------------------------------------------------
% Input
%   physcon.grav      ! Gravitational acceleration (m/s2)
%   physcon.tfrz      ! Freezing point of water (K)
%   physcon.visc0     ! Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
%   physcon.Dh0       ! Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
%   physcon.Dv0       ! Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
%   physcon.Dc0       ! Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)
%   atmos.patm        ! Atmospheric pressure (Pa)
%   atmos.rhomol      ! Molar density mol/m3)
%   atmos.wind        ! Wind speed (m/s)
%   atmos.tair        ! Air temperature (K)
%   leaf.dleaf        ! Leaf dimension (m)
%   flux.tleaf        ! Leaf temperature (K)
%
% Output
%   flux.gbh          ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
%   flux.gbv          ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
%   flux.gbc          ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
% -------------------------------------------------------------------------

% --- Adjust diffusivity for temperature and pressure

fac = 101325 / atmos.patm * (atmos.tair / physcon.tfrz)^1.81;

visc = physcon.visc0 * fac; % Kinematic viscosity (m2/s)
Dh = physcon.Dh0 * fac;     % Molecular diffusivity, heat (m2/s)
Dv = physcon.Dv0 * fac;     % Molecular diffusivity, H2O (m2/s)
Dc = physcon.Dc0 * fac;     % Molecular diffusivity, CO2 (m2/s)

% --- Dimensionless numbers

Re = atmos.wind * leaf.dleaf / visc; % Reynolds number
Pr = visc / Dh;                      % Prandtl number
Scv = visc / Dv;                     % Schmidt number for H2O
Scc = visc / Dc;                     % Schmidt number for CO2

% Grashof number

Gr = physcon.grav * leaf.dleaf^3 * max(flux.tleaf-atmos.tair, 0) / (atmos.tair * visc * visc);

% --- Empirical correction factor for Nu and Sh

b1 = 1.5;

% --- Nusselt number (Nu) and Sherwood numbers (H2O: Shv, CO2: Shc)

% Forced convection - laminar flow

Nu_lam  = b1 * 0.66 *  Pr^0.33 * Re^0.5;     % Nusselt number
Shv_lam = b1 * 0.66 * Scv^0.33 * Re^0.5;     % Sherwood number, H2O
Shc_lam = b1 * 0.66 * Scc^0.33 * Re^0.5;     % Sherwood number, CO2

% Forced convection - turbulent flow

Nu_turb  = b1 * 0.036 *  Pr^0.33 * Re^0.8;   % Nusselt number
Shv_turb = b1 * 0.036 * Scv^0.33 * Re^0.8;   % Sherwood number, H2O
Shc_turb = b1 * 0.036 * Scc^0.33 * Re^0.8;   % Sherwood number, CO2

% Choose correct flow regime for forced convection

Nu_forced = max(Nu_lam, Nu_turb);
Shv_forced = max(Shv_lam, Shv_turb);
Shc_forced = max(Shc_lam, Shc_turb);

% Free convection

Nu_free  = 0.54 *  Pr^0.25 * Gr^0.25;        % Nusselt number
Shv_free = 0.54 * Scv^0.25 * Gr^0.25;        % Sherwood number, H2O
Shc_free = 0.54 * Scc^0.25 * Gr^0.25;        % Sherwood number, CO2

% Both forced and free convection regimes occur together

Nu = Nu_forced + Nu_free;
Shv = Shv_forced + Shv_free;
Shc = Shc_forced + Shc_free;

% --- Boundary layer conductances (m/s)

flux.gbh = Dh *  Nu / leaf.dleaf;
flux.gbv = Dv * Shv / leaf.dleaf;
flux.gbc = Dc * Shc / leaf.dleaf;

% --- Convert conductance (m/s) to (mol/m2/s)

flux.gbh = flux.gbh * atmos.rhomol;
flux.gbv = flux.gbv * atmos.rhomol;
flux.gbc = flux.gbc * atmos.rhomol;
