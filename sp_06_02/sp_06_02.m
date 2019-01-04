% Supplemental program 6.2

% -----------------------------------------------------------
% Use Harman & Finnigan (2007, 2008) roughness sublayer (RSL)
% theory to obtain roughness lengths for momentum and scalars
% -----------------------------------------------------------

% --- Input parameters

% von Karman constant

vkc = 0.4;

% Leaf Nusselt number (heat) or Stanton number (scalar)
% with values of 0.1-0.2

%rc = 0.1;
rc = 0.2;

% Leaf drag coefficient

cd = 0.25;

% Canopy height (m)

hc = 20;

% Leaf area index (m2/m2)

LAI = 5;

% Leaf area density

lad = LAI / hc;

% Canopy density length scale (m)

Lc = 1 / (cd * lad);

% Obukhov length (m)

obu = -1000;

% --- Determine beta_val = u* / u(h) for the current Obukhov length

% Neutral value for beta = u* / u(h)

beta_neutral = 0.35;

% Lc/obu

LcL = Lc/obu;

% The unstable case is a quadratic equation for beta^2 at LcL

if (LcL <= 0)
   a = 1;
   b = 16 * LcL * beta_neutral^4;
   c = -beta_neutral^4;
   beta_val = sqrt((-b + sqrt(b^2 - 4 * a * c))/ (2 * a));

   % Error check

   y = beta_val^2 * LcL;
   fy = (1 - 16 * y)^(-0.25);
   err = beta_val * fy - beta_neutral;
   if (abs(err) > 1e-10)
      error('unstable case: error in beta')
   end
end

% The stable case is a cubic equation for beta at LcL

if (LcL > 0)
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

% --- For current beta = u*/u(h) determine displacement height

dp = beta_val^2 * Lc;                  % dp = hc - disp
disp = max(hc - dp, 0);                % Displacement height (m)

% Save canopy height (relative to displacement height),
% because this is used many time

h_minus_d = hc - disp;

% --- Turbulent Prandlt number (Pr) at canopy height

Prn = 0.5;         % Neutral value for Pr
Prvr = 0.3;        % Magnitude of variation of Pr with stability
Prsc = 2.0;        % Scale of variation of Pr with stability

Pr = Prn + Prvr * tanh(Prsc*Lc/obu);

% --- The "f" parameter relates the length scale of the scalar (heat) to that of momentum 

fval = (sqrt(1 + 4 * rc * Pr) - 1) / 2;

% --- Calculate the parameters c1 and c2 needed for the RSL function phi_hat

% Evaluate Monin-Obukhov phi functions at (hc-disp)/obu

[phi_m_hc] = phi_m_monin_obukhov (h_minus_d / obu);
[phi_c_hc] = phi_c_monin_obukhov (h_minus_d / obu);

% Roughness sublayer depth scale multiplier (dimensionless)

c2 = 0.5;

% c1 for momentum and scalars (dimensionless)

c1m = (1 -    vkc / (2 * beta_val * phi_m_hc)) * exp(c2/2);
c1c = (1 - Pr*vkc / (2 * beta_val * phi_c_hc)) * exp(c2/2);

% --- Evaluate the roughness sublayer psi_hat functions for momentum and scalars

% These are calculated at the canopy height. Note that here the heights are adjusted
% for the displacement height before the integration.

[psi_m_rsl_hc] = psi_m_rsl (h_minus_d, h_minus_d, obu, c1m, c2);  % momentum at (hc-disp)
[psi_c_rsl_hc] = psi_c_rsl (h_minus_d, h_minus_d, obu, c1c, c2);  % scalars at (hc-disp)

% --- Evaluate the Monin-Obukhov psi functions for momentum and scalars at the canopy height

[psi_m_hc] = psi_m_monin_obukhov (h_minus_d / obu);    % momentum at (hc-disp)/obu
[psi_c_hc] = psi_c_monin_obukhov (h_minus_d / obu);    % scalars at (hc-disp)/obu

% --- Roughness lengths z0m and z0c (m)

% z0m - Use bisection to find z0m, which lies between aval and bval, and refine the
% estimate until the difference is less than err

aval = hc;
bval = 0;
err = 1e-12;

[psi_m_z0m] = psi_m_monin_obukhov (aval / obu);
z0m = h_minus_d * exp(-vkc/beta_val) * exp(-psi_m_hc + psi_m_z0m) * exp(psi_m_rsl_hc);
fa = z0m - aval;

[psi_m_z0m] = psi_m_monin_obukhov (bval / obu);
z0m = h_minus_d * exp(-vkc/beta_val) * exp(-psi_m_hc + psi_m_z0m) * exp(psi_m_rsl_hc);
fb = z0m - bval;

if (fa * fb > 0)
   error('RSL bisection error: f(a) and f(b) do not have opposite signs')
end

while (abs(bval-aval) > err)
   cval = (aval + bval) / 2;
   [psi_m_z0m] = psi_m_monin_obukhov (cval / obu);
   z0m = h_minus_d * exp(-vkc/beta_val) * exp(-psi_m_hc + psi_m_z0m) * exp(psi_m_rsl_hc);
   fc = z0m - cval;
   if (fa * fc < 0)
      bval = cval; fb = fc;
   else
      aval = cval; fa = fc;
   end
end

z0m = cval;

% z0c - Use bisection to find z0c, which lies between aval and bval, and refine the
% estimate until the difference is less than err

aval = hc;
bval = 0;

[psi_c_z0c] = psi_c_monin_obukhov (aval / obu);
z0c = h_minus_d * exp(-vkc/beta_val*Pr/fval) * exp(-psi_c_hc + psi_c_z0c) * exp(psi_c_rsl_hc);
fa = z0c - aval;

[psi_c_z0c] = psi_c_monin_obukhov (bval / obu);
z0c = h_minus_d * exp(-vkc/beta_val*Pr/fval) * exp(-psi_c_hc + psi_c_z0c) * exp(psi_c_rsl_hc);
fb = z0c - bval;

if (fa * fb > 0)
   error('RSL bisection error: f(a) and f(b) do not have opposite signs')
end

while (abs(bval-aval) > err)
   cval = (aval + bval) / 2;
   [psi_c_z0c] = psi_c_monin_obukhov (cval / obu);
   z0c = h_minus_d * exp(-vkc/beta_val*Pr/fval) * exp(-psi_c_hc + psi_c_z0c) * exp(psi_c_rsl_hc);
   fc = z0c - cval;
   if (fa * fc < 0)
      bval = cval; fb = fc;
   else
      aval = cval; fa = fc;
   end
end

z0c = cval;

% --- Write output

kbinv = log(z0m/z0c);

fprintf('z0m = %15.3f\n',z0m)
fprintf('z0c = %15.3f\n',z0c)
fprintf('kB^{-1} = %15.3f\n',kbinv)
