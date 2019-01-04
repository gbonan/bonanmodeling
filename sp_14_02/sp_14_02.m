% Supplemental program 14.2

% ---------------------------------------------------------
% Test two-stream model with clumping. Compare numerical
% integration of leaf fluxes with canopy integrated fluxes.
% ---------------------------------------------------------

% ================
% Input parameters
% ================

% --- Leaf optical properties

rho_leaf = 0.10;                    % Leaf reflectance
tau_leaf = 0.05;                    % Leaf transmittance
omega_leaf = rho_leaf + tau_leaf;   % Leaf scattering coefficient
Kb = 0.58;                          % Direct beam extinction coefficient
Kd = 0.70;                          % Diffuse extinction coefficient
beta = 0.54;                        % Upscatter parameter for diffuse radiation
beta0 = 0.46;                       % Upscatter parameter for direct beam radiation

% --- Canopy variables

LAI = 6;                            % Leaf area index (m2/m2)
OMEGA = 0.75;                       % Clumping index

% --- Soil variables

albsoib = 0.1;                      % Soil albedo (direct beam)
albsoid = 0.1;                      % Soil albedo (diffuse)

% --- Atmospheric solar radiation, given here as a unit of incoming radiation

swskyb = 0.8;                       % Direct beam (W/m2)
swskyd = 0.2;                       % Diffuse (W/m2)

% ===================
% Analytical solution
% ===================

% --- Common terms: Eqs. (14.87) - (14.91)

b = (1 - (1 - beta) * omega_leaf) * Kd;
c = beta * omega_leaf * Kd;
h = sqrt(b*b - c*c);
u = (h - b - c) / (2 * h);
v = (h + b + c) / (2 * h);
g1 = (beta0 * Kb - b * beta0 - c * (1 - beta0)) * omega_leaf * Kb * swskyb / (h*h - Kb*Kb);
g2 = ((1 - beta0) * Kb + c * beta0 + b * (1 - beta0)) * omega_leaf * Kb * swskyb / (h*h - Kb*Kb);

% --- Exponential functions of leaf area

s1 = @(x) exp(-h * OMEGA * x);
s2 = @(x) exp(-Kb * OMEGA * x);

% --- Direct beam solution

% n1 (Eq. 14.92) and n2 (14.93)

num1 = v * (g1 + g2 * albsoid + albsoib * swskyb) * s2(LAI);
num2 = g2 * (u + v * albsoid) * s1(LAI);
den1 = v * (v + u * albsoid) / s1(LAI);
den2 = u * (u + v * albsoid) * s1(LAI);
n2b = (num1 - num2) / (den1 - den2);
n1b = (g2 - n2b * u) / v;

% Scattered direct beam fluxes:
% iupwb - direct beam flux scattered upward above cumulative LAI (W/m2); Eq. (14.94)
% idwnb - direct beam flux scattered downward below cumulative LAI (W/m2); Eq. (14.95)
% and their derivatives with respect to LAI

iupwb = @(x) -g1 * s2(x) + n1b * u * s1(x) + n2b * v / s1(x);
idwnb = @(x)  g2 * s2(x) - n1b * v * s1(x) - n2b * u / s1(x);
diupwb = @(x) ( Kb * g1 .* s2(x) - h * n1b * u .* s1(x) + h * n2b * v ./ s1(x)) * OMEGA;
didwnb = @(x) (-Kb * g2 .* s2(x) + h * n1b * v .* s1(x) - h * n2b * u ./ s1(x)) * OMEGA;

% icb - direct beam flux absorbed by canopy (W/m2); Eq. (14.97)

icb = swskyb * (1 - s2(LAI)) - iupwb(0) + iupwb(LAI) - idwnb(LAI);

% icsunb - direct beam flux absorbed by sunlit canopy (W/m2); Eq. (14.114)
% icshab - direct beam flux absorbed by shaded canopy (W/m2); Eq. (14.115)

a1b = -g1 *      (1 - s2(LAI)*s2(LAI)) / (2 * Kb) + ...
       n1b * u * (1 - s2(LAI)*s1(LAI)) / (Kb + h) + n2b * v * (1 - s2(LAI)/s1(LAI)) / (Kb - h);
a2b =  g2 *      (1 - s2(LAI)*s2(LAI)) / (2 * Kb) - ...
       n1b * v * (1 - s2(LAI)*s1(LAI)) / (Kb + h) - n2b * u * (1 - s2(LAI)/s1(LAI)) / (Kb - h);

icsunb = (1 - omega_leaf) * ((1 - s2(LAI)) * swskyb + Kd * (a1b + a2b) * OMEGA);
icshab = icb - icsunb;

% --- Diffuse solution

% n1 (Eq. 14.99) and n2 (14.100)

num = swskyd * (u + v * albsoid) * s1(LAI);
den1 = v * (v + u * albsoid) / s1(LAI);
den2 = u * (u + v * albsoid) * s1(LAI);
n2d = num / (den1 - den2);
n1d = -(swskyd + n2d * u) / v;

% Scattered diffuse fluxes:
% iupwd - diffuse flux scattered upward above cumulative LAI (W/m2); Eq. (14.101)
% idwnd - diffuse flux scattered downward below cumulative LAI (W/m2); Eq. (14.102)
% and their derivatives with respect to LAI

iupwd = @(x)  n1d * u * s1(x) + n2d * v / s1(x);
idwnd = @(x) -n1d * v * s1(x) - n2d * u / s1(x);
diupwd = @(x) (-h * n1d * u .* s1(x) + h * n2d * v ./ s1(x)) * OMEGA;
didwnd = @(x) ( h * n1d * v .* s1(x) - h * n2d * u ./ s1(x)) * OMEGA;

% icd - diffuse flux absorbed by canopy (W/m2); Eq. (14.104)

icd = swskyd - iupwd(0) + iupwd(LAI) - idwnd(LAI);

% icsund - diffuse flux absorbed by sunlit canopy (W/m2); Eq. (14.118)
% icshad - diffuse flux absorbed by shaded canopy (W/m2); Eq. (14.119)

a1d =  n1d * u * (1 - s2(LAI)*s1(LAI)) / (Kb + h) + n2d * v * (1 - s2(LAI)/s1(LAI)) / (Kb - h);
a2d = -n1d * v * (1 - s2(LAI)*s1(LAI)) / (Kb + h) - n2d * u * (1 - s2(LAI)/s1(LAI)) / (Kb - h);

icsund = (1 - omega_leaf) * Kd * (a1d + a2d) * OMEGA;
icshad = icd - icsund;

% --- Total canopy flux

ic = icb + icd;
icsun = icsunb + icsund;
icsha = icshab + icshad;

% ==================
% Numerical solution
% ==================

% --- Leaf fluxes for numerical solution

% fsun - sunlit fraction at cumulative LAI; Eq. (14.18)

fsun = @(x) OMEGA * exp(-Kb * OMEGA * x);

% ilbs - absorbed direct beam flux (scattered direct component) per leaf area
% at cumulative LAI, average for all leaves (J / m2 leaf / s); Eq. (14.108)

ilbs = @(x) omega_leaf * Kb * OMEGA * swskyb .* s2(x) + (diupwb(x) - didwnb(x));

% ild - absorbed diffuse flux per leaf area at cumulative LAI,
% average for all leaves (J / m2 leaf / s); Eq. (14.107)

ild = @(x) diupwd(x) - didwnd(x);

% ilsun - total absorbed flux (sunlit leaves) per sunlit leaf area
% at cumulative LAI (J / m2 leaf / s); Eq. (14.109)

ilsun = @(x) ild(x) + ilbs(x) + (1 - omega_leaf) * Kb * swskyb;

% ilsha - total absorbed flux (shaded leaves) per shaded leaf area
% at cumulative LAI (J / m2 leaf / s); Eq. (14.106)

ilsha = @(x) ild(x) + ilbs(x);

% il - total absorbed flux (average leaf) per unit leaf area
% at cumulative LAI (J / m2 leaf / s)

il = @(x) ilsun(x) .* fsun(x) + ilsha(x) .* (1 - fsun(x));

% --- Canopy fluxes (numerical)

% ic - absorbed solar radiation, total canopy (W/m2)
% icsun - absorbed solar radiation, sunlit canopy (W/m2); Eq. (14.110)
% icsha - absorbed solar radiation, shaded canopy (W/m2); Eq. (14.111)

eq_sun = @(x) ilsun(x) .* fsun(x);
eq_sha = @(x) ilsha(x) .* (1 - fsun(x));

numerical.ic = integral(il, 0, LAI);
numerical.icsun = integral(eq_sun, 0, LAI);
numerical.icsha = integral(eq_sha, 0, LAI);

% ============
% Print output
% ============

fprintf('ic = %15.5f %15.5f\n', ic, numerical.ic)
fprintf('icsun+icsha = %15.5f %15.5f\n', icsun+icsha, numerical.icsun+numerical.icsha)
fprintf('icsun = %15.5f %15.5f\n', icsun, numerical.icsun)
fprintf('icsha = %15.5f %15.5f\n', icsha, numerical.icsha)

% =====================================================================
% Check flux derivatives. Compare derivatives of upward and downward
% scattered fluxes with respect to leaf area used in the the two-stream
% equations (Eqs. 14.78, Eq. 14.79) with the derivatives obtained from
% the integrated analytical solution.
% =====================================================================

% Cumulative leaf area index to evaluate equations

x = 0.5;

% Direct beam: upward and downward fluxes

diupwb_two = (1 - (1 - beta) * omega_leaf) * Kd * OMEGA * iupwb(x) ...
           - beta * omega_leaf * Kd * OMEGA * idwnb(x) ...
           - beta0 * omega_leaf * Kb * OMEGA * swskyb * s2(x);
didwnb_two = -(1 - (1 - beta) * omega_leaf) * Kd * OMEGA * idwnb(x) ...
           + beta * omega_leaf * Kd * OMEGA * iupwb(x) ...
           + (1 - beta0) * omega_leaf * Kb * OMEGA * swskyb * s2(x);

fprintf('diupwb = %15.5f %15.5f\n', diupwb_two, diupwb(x))
fprintf('didwnb = %15.5f %15.5f\n', didwnb_two, didwnb(x))

% Diffuse: upward and downward fluxes

diupwd_two = (1 - (1 - beta) * omega_leaf) * Kd * OMEGA * iupwd(x) ...
           - beta * omega_leaf * Kd * OMEGA * idwnd(x);
didwnd_two = -(1 - (1 - beta) * omega_leaf) * Kd * OMEGA * idwnd(x) ...
           + beta * omega_leaf * Kd * OMEGA * iupwd(x);

fprintf('diupwd = %15.5f %15.5f\n', diupwd_two, diupwd(x))
fprintf('didwnd = %15.5f %15.5f\n', didwnd_two, didwnd(x))
