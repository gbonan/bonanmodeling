function [x] = leaf_boundary_layer (p)

% ----------------------------------------------------------------------------
% Calculate leaf boundary layer conductances
% ----------------------------------------------------------------------------

% Global variables

global tfrz g visc0 Dh0 Dv0 Dc0
global gbh gbw gbc dleaf
global pref tair wind rhomol 
global tleaf

% Adjust diffusivity for temperature and pressure

fac = 101325 / pref(p) * (tair(p) / tfrz)^1.81;

visc = visc0 * fac;    % Kinematic viscosity (m2/s)
Dh = Dh0 * fac;        % Molecular diffusivity, heat (m2/s)
Dv = Dv0 * fac;        % Molecular diffusivity, H2O (m2/s)
Dc = Dc0 * fac;        % Molecular diffusivity, CO2 (m2/s)

% Dimensionless numbers

Re = wind(p) * dleaf(p) / visc;   % Reynolds number
Pr  = visc / Dh;                  % Prandtl number
Scv = visc / Dv;                  % Schmidt number for H2O
Scc = visc / Dc;                  % Schmidt number for CO2
Gr = g * dleaf(p)^3 * max(tleaf(p) - tair(p), 0) / (tair(p) * visc * visc); % Grashof number

% Empirical correction factor for Nu and Sh

b1 = 1.5;

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

% Both convection regimes occur together

Nu = Nu_forced + Nu_free;
Shv = Shv_forced + Shv_free;
Shc = Shc_forced + Shc_free;

% Boundary layer conductances (mol/m2/s)

gbh(p) = Dh * Nu / dleaf(p) * rhomol(p);      % Heat
gbw(p) = Dv * Shv / dleaf(p) * rhomol(p);     % H2O
gbc(p) = Dc * Shc / dleaf(p) * rhomol(p);     % CO2

% Dummy output

x = 0;
