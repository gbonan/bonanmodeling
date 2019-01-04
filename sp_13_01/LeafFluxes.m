function [flux] = LeafFluxes (physcon, atmos, leaf, flux)

% --- Leaf temperature, energy fluxes, photosynthesis, and stomatal conductance
% using stomatal optimization

[flux] = StomataOptimization (physcon, atmos, leaf, flux);
