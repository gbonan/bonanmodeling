function [fluxvar, fx] = rsl (physcon, forcvar, surfvar, fluxvar, x)

% Use Harman & Finnigan (2007, 2008) roughness sublayer (RSL) theory to obtain
% the Obukhov length (obu). This is the function to solve for the Obukhov
% length. For the current estimate of the Obukhov length (x), calculate
% u*, T*, and q* and then the new length (obu). The function value is the
% change in Obukhov length: fx = x - obu.

% -------------------------------------------------------------------------
% Input
%   x                   ! Current estimate for Obukhov length (m)
%   physcon.vkc         ! von Karman constant
%   physcon.grav        ! Gravitational acceleration (m/s2)
%   physcon.mmh2o       ! Molecular mass of water (kg/mol)
%   forcvar.zref        ! Reference height (m)
%   forcvar.uref        ! Wind speed at reference height (m/s)
%   forcvar.thref       ! Potential temperature at reference height (K)
%   forcvar.thvref      ! Virtual potential temperature at reference height (K)
%   forcvar.eref        ! Vapor pressure at reference height (Pa)
%   forcvar.pref        ! Atmospheric pressure (Pa)
%   forcvar.mmair       ! Molecular mass of air at reference height (kg/mol)
%   fluxvar.tsrf        ! Surface temperature (K)
%   fluxvar.esrf        ! Surface vapor pressure (Pa)
%   surfvar.Lc          ! Canopy density length scale (m)
%   surfvar.hc          ! Canopy height (m)
%   surfvar.rc          ! Leaf Nusselt number (heat) or Stanton number (scalar)
%
% Output
%   fluxvar.z0m         ! Roughness length for momentum (m)
%   fluxvar.z0c         ! Roughness length for scalars (m)
%   fluxvar.disp        ! Displacement height (m)
%   fluxvar.ustar       ! Friction velocity (m/s)
%   fluxvar.tstar       ! Temperature scale (K)
%   fluxvar.qstar       ! Water vapor scale (mol/mol)
%   fluxvar.obu         ! Obukhov length (m)
%   fx                  ! Change in Obukhov length (x - obu)
% -------------------------------------------------------------------------

% --- Prevent near-zero values of Obukhov length

if (abs(x) <= 0.1)
   x = 0.1;
end

% --- Determine beta_val = u* / u(h) for the current Obukhov length

% Neutral value for beta = u* / u(h)

beta_neutral = 0.35;

% Lc/obu

LcL = surfvar.Lc/x;

if (LcL <= 0)

   % The unstable case is a quadratic equation for beta^2 at LcL

   a = 1;
   b = 16 * LcL * beta_neutral^4;
   c = -beta_neutral^4;
   beta_val = sqrt((-b + sqrt(b^2 - 4 * a * c)) / (2 * a));

   % Error check

   y = beta_val^2 * LcL;
   fy = (1 - 16 * y)^(-0.25);
   err = beta_val * fy - beta_neutral;
   if (abs(err) > 1e-10)
      error('unstable case: error in beta')
   end

else

   % The stable case is a cubic equation for beta at LcL

   a = 5 * LcL;
   b = 0;
   c = 1;
   d = -beta_neutral;
   q = (2*b^3 - 9*a*b*c + 27*(a^2)*d)^2 - 4*(b^2 - 3*a*c)^3;
   q = sqrt(q);
   r = 0.5 * (q + 2*b^3 - 9*a*b*c + 27*(a^2)*d);
   r = r^(1/3);
   beta_val = -(b+r)/(3*a) - (b^2 - 3*a*c)/(3*a*r);

   % Error check

   y = beta_val^2 * LcL;
   fy = 1 + 5 * y;
   err = beta_val * fy - beta_neutral;
   if (abs(err) > 1e-10)
      error('stable case: error in beta')
   end

end

% Place limits on beta

beta_val = min(0.5, max(beta_val,0.2));

% --- For current beta = u*/u(h) determine displacement height

dp = beta_val^2 * surfvar.Lc;                  % dp = hc - disp
fluxvar.disp = max(surfvar.hc - dp, 0);        % Displacement height (m)

% Save reference height and canopy height (relative to displacement height),
% because these are used many times

z_minus_d = forcvar.zref - fluxvar.disp;
h_minus_d = surfvar.hc - fluxvar.disp;

% --- Turbulent Prandlt number (Pr) at canopy height

Prn = 0.5;         % Neutral value for Pr
Prvr = 0.3;        % Magnitude of variation of Pr with stability
Prsc = 2.0;        % Scale of variation of Pr with stability

Pr = Prn + Prvr * tanh(Prsc*surfvar.Lc/x);

% --- The "f" parameter relates the length scale of the scalar (heat) to that of momentum 

f = (sqrt(1 + 4 * surfvar.rc * Pr) - 1) / 2;

% --- Calculate the parameters c1 and c2 needed for the RSL function phi_hat

% Evaluate Monin-Obukhov phi functions at (hc-disp)/obu

[phi_m_hc] = phi_m_monin_obukhov (h_minus_d / x);
[phi_c_hc] = phi_c_monin_obukhov (h_minus_d / x);

% Roughness sublayer depth scale multiplier (dimensionless)

c2 = 0.5;

% c1 for momentum and scalars (dimensionless)

c1m = (1 -    physcon.vkc / (2 * beta_val * phi_m_hc)) * exp(c2/2);
c1c = (1 - Pr*physcon.vkc / (2 * beta_val * phi_c_hc)) * exp(c2/2);

% --- Evaluate the roughness sublayer psi_hat functions for momentum and scalars

% These are calculated at the reference height and at the canopy height. Note that
% here the heights are adjusted for the displacement height before the integration.

[psi_m_rsl_zref] = psi_m_rsl (z_minus_d, h_minus_d, x, c1m, c2);  % momentum at (zref-disp)
[psi_m_rsl_hc]   = psi_m_rsl (h_minus_d, h_minus_d, x, c1m, c2);  % momentum at (hc-disp)

[psi_c_rsl_zref] = psi_c_rsl (z_minus_d, h_minus_d, x, c1c, c2);  % scalars at (zref-disp)
[psi_c_rsl_hc]   = psi_c_rsl (h_minus_d, h_minus_d, x, c1c, c2);  % scalars at (hc-disp)

% --- Evaluate the Monin-Obukhov psi functions for momentum and scalars

% These are calculated at the reference height and at the canopy height

[psi_m_zref] = psi_m_monin_obukhov (z_minus_d / x);    % momentum at (zref-disp)/obu
[psi_m_hc]   = psi_m_monin_obukhov (h_minus_d / x);    % momentum at (hc-disp)/obu

[psi_c_zref] = psi_c_monin_obukhov (z_minus_d / x);    % scalars at (zref-disp)/obu
[psi_c_hc]   = psi_c_monin_obukhov (h_minus_d / x);    % scalars at (hc-disp)/obu

% --- Calculate u* (m/s), T* (K), q* (mol/mol), and Tv* (K)

zlog = log(z_minus_d / h_minus_d);
psim = -psi_m_zref + psi_m_hc + psi_m_rsl_zref - psi_m_rsl_hc + physcon.vkc / beta_val;
psic = -psi_c_zref + psi_c_hc + psi_c_rsl_zref - psi_c_rsl_hc + physcon.vkc / beta_val * Pr / f;

fluxvar.ustar = forcvar.uref * physcon.vkc / (zlog + psim);
fluxvar.tstar = (forcvar.thref - fluxvar.tsrf) * physcon.vkc / (zlog + psic);
fluxvar.qstar = (forcvar.eref - fluxvar.esrf) / forcvar.pref * physcon.vkc / (zlog + psic);
tvstar = fluxvar.tstar + 0.61 * forcvar.thref * fluxvar.qstar * (physcon.mmh2o / forcvar.mmair);

% --- Calculate Obukhov length (m)

fluxvar.obu = fluxvar.ustar^2 * forcvar.thvref / (physcon.vkc * physcon.grav * tvstar);
fx = x - fluxvar.obu;

% --- Roughness lengths z0m and z0c (m)

% z0m - Use bisection to find z0m, which lies between aval and bval, and refine the
% estimate until the difference is less than err

aval = surfvar.hc;
bval = 0;
err = 1e-08;

[psi_m_z0m] = psi_m_monin_obukhov (aval / x);
z0m = h_minus_d * exp(-physcon.vkc/beta_val) * exp(-psi_m_hc + psi_m_z0m) * exp(psi_m_rsl_hc);
fa = z0m - aval;

[psi_m_z0m] = psi_m_monin_obukhov (bval / x);
z0m = h_minus_d * exp(-physcon.vkc/beta_val) * exp(-psi_m_hc + psi_m_z0m) * exp(psi_m_rsl_hc);
fb = z0m - bval;

if (fa * fb > 0)
   error('RSL bisection error: f(a) and f(b) do not have opposite signs')
end

while (abs(bval-aval) > err)
   cval = (aval + bval) / 2;

   [psi_m_z0m] = psi_m_monin_obukhov (cval / x);
   z0m = h_minus_d * exp(-physcon.vkc/beta_val) * exp(-psi_m_hc + psi_m_z0m) * exp(psi_m_rsl_hc);
   fc = z0m - cval;

   if (fa * fc < 0)
      bval = cval; fb = fc;
   else
      aval = cval; fa = fc;
   end
end

fluxvar.z0m = cval;

% z0c - Use bisection to find z0c, which lies between aval and bval, and refine the
% estimate until the difference is less than err

aval = surfvar.hc;
bval = 0;

[psi_c_z0c] = psi_c_monin_obukhov (aval / x);
z0c = h_minus_d * exp(-physcon.vkc/beta_val*Pr/f) * exp(-psi_c_hc + psi_c_z0c) * exp(psi_c_rsl_hc);
fa = z0c - aval;

[psi_c_z0c] = psi_c_monin_obukhov (bval / x);
z0c = h_minus_d * exp(-physcon.vkc/beta_val*Pr/f) * exp(-psi_c_hc + psi_c_z0c) * exp(psi_c_rsl_hc);
fb = z0c - bval;

if (fa * fb > 0)
   error('RSL bisection error: f(a) and f(b) do not have opposite signs')
end

while (abs(bval-aval) > err)
   cval = (aval + bval) / 2;

   [psi_c_z0c] = psi_c_monin_obukhov (cval / x);
   z0c = h_minus_d * exp(-physcon.vkc/beta_val*Pr/f) * exp(-psi_c_hc + psi_c_z0c) * exp(psi_c_rsl_hc);
   fc = z0c - cval;

   if (fa * fc < 0)
      bval = cval; fb = fc;
   else
      aval = cval; fa = fc;
   end
end

fluxvar.z0c = cval;
