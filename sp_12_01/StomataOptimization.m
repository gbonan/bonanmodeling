function [flux] = StomataOptimization (physcon, atmos, leaf, flux)

% Leaf temperature, energy fluxes, photosynthesis, and stomatal conductance
% with water-use efficiency stomatal optimization

% --- Low and high initial estimates for gs (mol H2O/m2/s)

gs1 = 0.002;
gs2 = 2.0;

% --- Check for minimum stomatal conductance linked to low light
% based on the efficiency check for gs1 and gs2 (check1, check2)

[flux, check1] = StomataEfficiency (physcon, atmos, leaf, flux, gs1);
[flux, check2] = StomataEfficiency (physcon, atmos, leaf, flux, gs2);

if (check1 * check2 < 0)

   % Calculate gs using the function StomataEfficiency to iterate gs
   % to an accuracy of tol (mol H2O/m2/s)

   tol = 0.004;
   func_name = 'StomataEfficiency';
   [flux, root] = brent_root (func_name, physcon, atmos, leaf, flux, gs1, gs2, tol);
   flux.gs = root;

else

   % Low light. Set gs to minimum conductance

   flux.gs = 0.002;

end

% --- Leaf fluxes for this gs

[flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);
[flux] = LeafTranspiration (physcon, atmos, flux);
