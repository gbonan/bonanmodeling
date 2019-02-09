function [psi_hat_m] = psi_m_rsl (z, h, L, c1, c2)

% --- Evaluate the roughness sublayer (RSL) function psi_hat for momentum
% at z. Note that z has already been adjusted for the displacement height
% (i.e., using z - d).

% ------------------------------------------------------
% Input
%   z            ! Vertical height - displacement height (m)
%   h            ! Canopy height - displacement height (m)
%   L            ! Obukhov length (m)
%   c1           ! Parameter for RSL function phi_hat (dimensionless)
%   c2           ! Parameter for RSL function phi_hat (dimensionless)
%
% Output
%   psi_hat_m    ! RSL psi_hat function for momentum (dimensionless)
% ------------------------------------------------------

% The function to integrate depends on unstable (f1) or stable (f2)

f1 = @(x) (1-16*x/L).^(-0.25) .* (1 - (1 - c1*exp(-c2*x/(2*h)))) ./ x;
f2 = @(x) (1+5*x/L)           .* (1 - (1 - c1*exp(-c2*x/(2*h)))) ./ x;

% Numerically integrate the function from z to infinity

if (L < 0)
   psi_hat_m = integral (f1, z, inf);
else
   psi_hat_m = integral (f2, z, inf);
end
