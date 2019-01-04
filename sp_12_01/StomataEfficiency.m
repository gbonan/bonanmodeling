function [flux, val] = StomataEfficiency (physcon, atmos, leaf, flux, gs_val)

% Stomatal efficiency check to determine maximum gs. For the stomatal conductance
% gs_val, calculate photosynthesis for an increase in stomatal conductance equal
% to "delta". The returned value is positive if this increase produces a change
% in photosynthesis > iota*vpd*delta. The returned value is negative if the increase
% produces a change in photosynthesis < iota*vpd*delta.

% --- Specify "delta" as a small difference in gs (mol H2O/m2/s)

delta = 0.001;

% -- Calculate photosynthesis at lower gs (gs_val - delta)
% gs2 - lower value for gs (mol H2O/m2/s)
% an2 - leaf photosynthesis at gs2 (umol CO2/m2/s)

gs2 = gs_val - delta;
flux.gs = gs2;
[flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);
an2 = flux.an;

% -- Calculate photosynthesis at higher gs (gs_val)
% gs1 - higher value for gs (mol H2O/m2/s)
% an1 - leaf photosynthesis at gs1 (umol CO2/m2/s)

gs1 = gs_val;
flux.gs = gs1;
[flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);
an1 = flux.an;

% --- Efficiency check: val < 0 when d(An) / d(gs) < iota * vpd

val = (an1 - an2) - leaf.iota * delta * (flux.vpd / atmos.patm);
