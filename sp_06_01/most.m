function [fx] = most (x, var)

% -------------------------------------------------------------------------
% Use Monin-Obukhov similarity theory to obtain the Obukhov length (L).
%
% This is the function to solve for the Obukhov length. For current estimate
% of the Obukhov length (x), calculate u* and T* and then the new length (L).
% The function value is the change in Obukhov length: fx = x - L.
%
% Input:  x        ! Current estimate for Obukhov length (m)
%         var.z1   ! Height (m)
%         var.z2   ! Height (m)
%         var.u1   ! Wind speed at z1 (m/s)
%         var.u2   ! Wind speed at z2 (m/s)
%         var.t1   ! Temperature at z1 (m/s)
%         var.t2   ! Temperature at z2 (m/s)
%         var.d    ! Displacement height (m)
%         var.k    ! von Karman constant
%         var.g    ! Gravitational acceleration (m/s2)
% Output: fx       ! Change in Obukhov length (x - L)
%
% Local:  psi_m_z2 ! psi for momentum at height z2 (dimensionless)
%         psi_m_z1 ! psi for momentum at height z1 (dimensionless)
%         psi_c_z2 ! psi for scalars at height z2 (dimensionless)
%         psi_c_z1 ! psi for scalars at height z1 (dimensionless)
%         ustar    ! Friction velocity (m/s)
%         tstar    ! Temperature scale (K)
%         L        ! Obukhov length (m)
% -------------------------------------------------------------------------

% Prevent near-zero values of Obukhov length

if (abs(x) <= 0.1)
   x = 0.1;
end

% Evaluate psi for momentum at heights z2 and z1

[psi_m_z2] = psi_m_monin_obukhov((var.z2-var.d)/x);
[psi_m_z1] = psi_m_monin_obukhov((var.z1-var.d)/x);

% Evaluate psi for scalars at heights z2 and z1

[psi_c_z2] = psi_c_monin_obukhov((var.z2-var.d)/x);
[psi_c_z1] = psi_c_monin_obukhov((var.z1-var.d)/x);

% Calculate u* (m/s) and T* (K)

ustar = (var.u2 - var.u1) * var.k / (log((var.z2-var.d)/(var.z1-var.d)) - (psi_m_z2 - psi_m_z1));
tstar = (var.t2 - var.t1) * var.k / (log((var.z2-var.d)/(var.z1-var.d)) - (psi_c_z2 - psi_c_z1));

% Calculate L (m)

L = ustar^2 * var.t2 / (var.k * var.g * tstar);

% Calculate change in L

fx = x - L;
