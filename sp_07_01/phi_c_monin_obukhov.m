function [phi_c] = phi_c_monin_obukhov (x)

% --- Evaluate the Monin-Obukhov phi function for scalars at x

if (x < 0)
   phi_c = (1 - 16 * x)^(-0.5);
else
   phi_c = 1 + 5 * x;
end
