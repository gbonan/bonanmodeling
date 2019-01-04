function [dtleaf] = leaf_temperature (p)

% ----------------------------------------------------------------------------
% Use Newton-Raphson iteration to solve the leaf energy budget for leaf temperature
% ----------------------------------------------------------------------------

% Global variables

global tfrz mmh2o sigma
global tair cpair pref eair
global gsw gbw gbh emleaf
global tleaf qa rn lwrad sh lh

niter = 0;           % Number of iterations
err = 1e36;          % Energy inbalance (W/m2)

% Iteration is until energy imbalance < 1e-06 W/m2 or to 100 iterations

while (niter <= 100 & abs(err) > 1e-06)

   % Increment iteration counter

   niter = niter + 1;

   % Saturation vapor pressure ESAT (Pa) and temperature derivative DESAT (Pa/K)

   [esat, desat] = satvap (tleaf(p)-tfrz);

   % Latent heat of vaporization (J/mol)

   lambda = 2501.6 - 2.3773 * (tair(p) - tfrz);  % J/g
   lambda = lambda * 1000 * mmh2o;               % J/g -> J/kg -> J/mol

   % Leaf conductance for water vapor (mol H2O/m2/s)

   gleaf = gsw(p) * gbw(p) / (gsw(p) + gbw(p));

   % Emitted longwave radiation LWRAD (W/m2) and temperature derivative DLWRAD (W/m2/K)

   lwrad(p) = 2 * emleaf(p) * sigma * tleaf(p)^4;
   dlwrad = -8 * emleaf(p) * sigma * tleaf(p)^3;

   % Sensible heat flux SH (W/m2) and temperature derivative DSH (W/m2/K)

   sh(p) = 2 * cpair(p) * (tleaf(p) - tair(p)) * gbh(p);
   dsh = -2 * cpair(p) * gbh(p);

   % Latent heat flux LH (W/m2) and temperature derivative DLH (W/m2/K)

   lh(p) = lambda * (esat - eair(p)) / pref(p) * gleaf;
   dlh = -lambda * desat / pref(p) * gleaf;

   % Energy imbalance (W/m2)

   err = qa(p) - lwrad(p) - sh(p) - lh(p);

   % Change in leaf temperature (K)

   dtleaf = -err / (dlwrad + dsh + dlh);

   % Update leaf temperature (K)

   tleaf(p) = tleaf(p) + dtleaf;

end

% Net radiation (W/m2)

rn(p) = qa(p) - lwrad(p);

% Error check

err = rn(p) - sh(p) - lh(p);
if (abs(err) > 1e-06)
   fprintf('err = %15.3f\n',err)
   fprintf('qa = %15.3f\n',qa(p))
   fprintf('lwrad = %15.3f\n',lwrad(p))
   fprintf('sh = %15.3f\n',sh(p))
   fprintf('lh = %15.3f\n',lh(p))
   error ('leaf temperature error')
end
