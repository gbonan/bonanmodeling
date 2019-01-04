function [leafwp] = LeafWaterPotential (physcon, leaf, flux)

% Calculate leaf water potential for the current transpiration rate

% -------------------------------------------------------------------------
% Input
%   physcon.grav     ! Gravitational acceleration (m/s2)
%   physcon.denh2o   ! Density of liquid water (kg/m3)
%   leaf.capac       ! Plant capacitance (mmol H2O/m2 leaf area/MPa)
%   flux.height      ! Leaf height (m)
%   flux.lsc         ! Leaf-specific conductance (mmol H2O/m2 leaf/s/MPa)
%   flux.psi_leaf    ! Leaf water potential at beginning of time step (MPa)
%   flux.psi_soil    ! Weighted soil water potential (MPa)
%   flux.etflx       ! Leaf transpiration flux (mol H2O/m2 leaf/s)
%   flux.dt          ! Model time step (s)
% Output
%   leafwp           ! Leaf water potential at end of time step (MPa)
% -------------------------------------------------------------------------

% --- Head of pressure (MPa/m)

head = physcon.denh2o * physcon.grav * 1e-06;

% --- Change in leaf water potential is: dy / dt = (a - y) / b. The integrated change 
% over a full model timestep is: dy = (a - y0) * (1 - exp(-dt/b))

y0 = flux.psi_leaf;
a = flux.psi_soil - head * flux.height - 1000 * flux.etflx / flux.lsc;
b = leaf.capac / flux.lsc;
dy = (a - y0) * (1 - exp(-flux.dt/b));
leafwp = y0 + dy;
