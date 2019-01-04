function [psi_hat_c] = psi_c_rsl (z, h, L, c1, c2)

% --- Evaluate the roughness sublayer (RSL) function psi_hat for scalars
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
%   psi_hat_c    ! RSL psi_hat function for scalars (dimensionless)
% ------------------------------------------------------

% The function to integrate depends on unstable (f1) or stable (f2)

f1 = @(x) (1-16*x/L).^(-0.5) .* (1 - (1 - c1*exp(-c2*x/(2*h)))) ./ x;
f2 = @(x) (1+5*x/L)          .* (1 - (1 - c1*exp(-c2*x/(2*h)))) ./ x;

% Numerically integrate the function from z to infinity

if (L < 0)
   psi_hat_c = integral (f1, z, inf);
else
   psi_hat_c = integral (f2, z, inf);
end
