function [flux, val] = StomataEfficiency (physcon, atmos, leaf, flux, gs_val)

% Stomata water-use efficiency check and cavitation check to determine maximum gs. 
% For the stomatal conductance gs_val, calculate photosynthesis and leaf
% water potential for an increase in stomatal conductance equal to "delta".
% The returned value is positive if this increase produces a change in
% photosynthesis > iota*vpd*delta or if the leaf water potential is > minlwp.
% The returned value is negative if the increase produces a change in
% photosynthesis < iota*vpd*delta or if the leaf water potential is < minlwp.

% --- Leaf boundary layer conductances

[flux] = LeafBoundaryLayer (physcon, atmos, leaf, flux);

% --- Specify "delta" as a small difference in gs (mol H2O/m2/s)

delta = 0.001;

% --- Calculate photosynthesis at lower gs (gs_val - delta), but first need leaf temperature for this gs
% gs2 - lower value for gs (mol H2O/m2/s)
% an2 - leaf photosynthesis at gs2 (umol CO2/m2/s)

gs2 = gs_val - delta;
flux.gs = gs2;
[flux] = LeafTemperature (physcon, atmos, leaf, flux);
[flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);
an2 = flux.an;

% --- Calculate photosynthesis at higher gs (gs_val), but first need leaf temperature for this gs
% gs1 - higher value for gs (mol H2O/m2/s)
% an1 - leaf photosynthesis at gs1 (umol CO2/m2/s)

gs1 = gs_val;
flux.gs = gs1;
[flux] = LeafTemperature (physcon, atmos, leaf, flux);
[flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);
an1 = flux.an;

% --- Efficiency check: wue < 0 when d(An) / d(gs) < iota * vpd

wue = (an1 - an2) - leaf.iota * delta * (flux.vpd / atmos.patm);

% --- Cavitation check: minpsi < 0 when leafwp < minlwp

[leafwp] = LeafWaterPotential (physcon, leaf, flux);
minpsi = leafwp - leaf.minlwp;

% --- Return the minimum of the two checks

val = min(wue, minpsi);
