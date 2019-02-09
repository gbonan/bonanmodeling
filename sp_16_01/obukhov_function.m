function [fluxvar, fx] = obukhov_function (physcon, forcvar, surfvar, fluxvar, x)

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
%   forcvar.qref        ! Water vapor at reference height (mol/mol)
%   forcvar.rhomol      ! Molar density (mol/m3)
%   forcvar.mmair       ! Molecular mass of air at reference height (kg/mol)
%   fluxvar.taf         ! Air temperature at canopy top (K)
%   fluxvar.qaf         ! Water vapor at canopy top (mol/mol)
%   fluxvar.Lc          ! Canopy density length scale (m)
%   surfvar.hc          ! Canopy height (m)
%   surfvar.p           ! Index for grid point to process
%
% Output
%   fluxvar.c1m         ! Roughness sublayer c1 parameter for momentum (dimensionless)
%   fluxvar.c1c         ! Roughness sublayer c1 parameter for scalars (dimensionless)
%   fluxvar.c2          ! Roughness sublayer depth scale multiplier (dimensionless)
%   fluxvar.disp        ! Displacement height (m)
%   fluxvar.beta        ! u* / u(hc)
%   fluxvar.PrSc        ! Prandtl (Schmidt) number at canopy top
%   fluxvar.ustar       ! Friction velocity (m/s)
%   fluxvar.tstar       ! Temperature scale (K)
%   fluxvar.qstar       ! Water vapor scale (mol/mol)
%   fluxvar.gac         ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
%   fluxvar.obu_ustar   ! Obukhov length used for u* (m)
%   fluxvar.obu         ! Value for Obukhov length (m)
%   fx                  ! Change in Obukhov length (x - obu)
% -------------------------------------------------------------------------

% --- Index for grid point to process

p = surfvar.p;

% --- Prevent near-zero values of Obukhov length

if (abs(x) <= 0.1)
   x = 0.1;
end

% --- Determine beta_val = u* / u(hc) for the current Obukhov length

% Neutral value for beta = u* / u(hc)

beta_neutral = 0.35;

% Lc/obu

LcL = fluxvar.Lc(p)/x;

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
      error('obukhov_function: unstable case - error in beta')
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
      error('obukhov_function: stable case - error in beta')
   end

end

% Place limits on beta

beta_val = min(0.5, max(beta_val,0.2));
fluxvar.beta(p) = beta_val;

% --- For current beta = u*/u(hc) determine displacement height

dp = beta_val^2 * fluxvar.Lc(p);               % dp = hc - disp
fluxvar.disp(p) = max(surfvar.hc(p) - dp, 0);  % Displacement height (m)

% Save reference height and canopy height (relative to displacement height)
% because these are used many times

z_minus_d = forcvar.zref(p) - fluxvar.disp(p);
h_minus_d = surfvar.hc(p) - fluxvar.disp(p);

% --- Turbulent Prandlt (Schmidt) number (PrSc) at canopy height

PrSc0 = 0.5;        % Neutral value for Pr (Sc)
PrSc1 = 0.3;        % Magnitude of variation of Pr (Sc) with stability
PrSc2 = 2.0;        % Scale of variation of Pr (Sc) with stability

fluxvar.PrSc(p) = PrSc0 + PrSc1 * tanh(PrSc2*fluxvar.Lc(p)/x);

% --- Calculate the parameters c1 and c2 needed for the RSL function phi_hat

% Evaluate Monin-Obukhov phi functions at (hc-disp)/obu

[phi_m_hc] = phi_m_monin_obukhov (h_minus_d / x);
[phi_c_hc] = phi_c_monin_obukhov (h_minus_d / x);

% Roughness sublayer depth scale multiplier (dimensionless)

fluxvar.c2 = 0.5;

% c1 for momentum and scalars (dimensionless)

fluxvar.c1m(p) = (1 -                 physcon.vkc / (2 * beta_val * phi_m_hc)) * exp(fluxvar.c2/2);
fluxvar.c1c(p) = (1 - fluxvar.PrSc(p)*physcon.vkc / (2 * beta_val * phi_c_hc)) * exp(fluxvar.c2/2);

% --- Evaluate the roughness sublayer psi_hat functions for momentum and scalars

% These are calculated at the reference height and at the canopy height. Note that
% here the heights are adjusted for the displacement height before the integration.

[psi_m_rsl_zref] = psi_m_rsl (z_minus_d, h_minus_d, x, fluxvar.c1m(p), fluxvar.c2);  % momentum at (zref-disp)
[psi_m_rsl_hc]   = psi_m_rsl (h_minus_d, h_minus_d, x, fluxvar.c1m(p), fluxvar.c2);  % momentum at (hc-disp)

[psi_c_rsl_zref] = psi_c_rsl (z_minus_d, h_minus_d, x, fluxvar.c1c(p), fluxvar.c2);  % scalars at (zref-disp)
[psi_c_rsl_hc]   = psi_c_rsl (h_minus_d, h_minus_d, x, fluxvar.c1c(p), fluxvar.c2);  % scalars at (hc-disp)

% --- Evaluate the Monin-Obukhov psi functions for momentum and scalars

% These are calculated at the reference height and at the canopy height

[psi_m_zref] = psi_m_monin_obukhov (z_minus_d / x);    % momentum at (zref-disp)/obu
[psi_m_hc]   = psi_m_monin_obukhov (h_minus_d / x);    % momentum at (hc-disp)/obu

[psi_c_zref] = psi_c_monin_obukhov (z_minus_d / x);    % scalars at (zref-disp)/obu
[psi_c_hc]   = psi_c_monin_obukhov (h_minus_d / x);    % scalars at (hc-disp)/obu

% --- Calculate u* (m/s), T* (K), q* (mol/mol), and Tv* (K)

zlog = log(z_minus_d / h_minus_d);
psim = -psi_m_zref + psi_m_hc + psi_m_rsl_zref - psi_m_rsl_hc + physcon.vkc / beta_val;
psic = -psi_c_zref + psi_c_hc + psi_c_rsl_zref - psi_c_rsl_hc;

fluxvar.ustar(p) = forcvar.uref(p) * physcon.vkc / (zlog + psim);
fluxvar.tstar(p) = (forcvar.thref(p) - fluxvar.taf(p)) * physcon.vkc / (zlog + psic);
fluxvar.qstar(p) = (forcvar.qref(p) - fluxvar.qaf(p)) * physcon.vkc / (zlog + psic);
fluxvar.obu_ustar(p) = x;

% --- Aerodynamic conductance to canopy height

fluxvar.gac(p) = forcvar.rhomol(p) * physcon.vkc * fluxvar.ustar(p) / (zlog + psic);

% --- Calculate new Obukhov length (m)

tvstar = fluxvar.tstar(p) + 0.61 * forcvar.thref(p) * fluxvar.qstar(p) * (physcon.mmh2o / forcvar.mmair(p));
fluxvar.obu(p) = fluxvar.ustar(p)^2 * forcvar.thvref(p) / (physcon.vkc * physcon.grav * tvstar);
fx = x - fluxvar.obu(p);
