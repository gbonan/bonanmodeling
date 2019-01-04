function [psi_c] = psi_c_monin_obukhov (x)

% --- Evaluate the Monin-Obukhov psi function for scalars at x

if (x < 0)
   y = (1 - 16 * x)^0.25;
   psi_c = 2 * log((1 + y^2)/2);
else
   psi_c = -5 * x;
end
