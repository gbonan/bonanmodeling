function [phi_m] = phi_m_monin_obukhov (x)

% --- Evaluate the Monin-Obukhov phi function for momentum at x

if (x < 0)
   phi_m = (1 - 16 * x)^(-0.25);
else
   phi_m = 1 + 5 * x;
end
