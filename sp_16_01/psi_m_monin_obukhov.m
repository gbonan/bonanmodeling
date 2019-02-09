function [psi_m] = psi_m_monin_obukhov (x)

% --- Evaluate the Monin-Obukhov psi function for momentum at x

if (x < 0)
   y = (1 - 16 * x)^0.25;
   psi_m = 2 * log((1 + y)/2) + log((1 + y^2)/2) - 2 * atan(y) + pi / 2;
else
   psi_m = -5 * x;
end
