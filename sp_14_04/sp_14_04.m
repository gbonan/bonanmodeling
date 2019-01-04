% Supplemental program 14.4

% -------------------------------------------------------------------
% Compare longwave radiative transfer calculated using the analytical
% canopy-integrated model and the Norman multilayer model
% -------------------------------------------------------------------

% --- Parameters

sigma = 5.67e-08;             % Stefan-Boltzmann constant (W/m2/K4)
tfrz = 273.15;                % Freezing point of water (K)
emleaf = 0.98;                % Leaf emissivity
emgrnd = 1.00;                % Ground (soil) emissivity
% emgrnd = 0.96;                % Ground (soil) emissivity

LAI = 4.9;                    % Leaf area index (m2/m2)
tveg = tfrz + 25;             % Canopy temperature (K)
tgrnd = tfrz + 20;            % Ground temperature (K)
irsky = 400;                  % Atmospheric longwave radiation (W/m2)

% --- For Norman radiation

nveg = 49;                    % Number of leaf layers (each with lai = dlai)
nsoi = 1;                     % First canopy layer is soil
nbot = nsoi + 1;              % Index for bottom leaf layer
ntop = nbot + nveg - 1;       % Index for top leaf layer

for iv = nbot:ntop
   tleaf(iv) = tveg;          % Leaf temperature (K)
   dlai(iv) = 0.1;            % Layer leaf area index (m2/m2)
   td(iv) = 0.915;            % Exponential transmittance of diffuse radiation through a single leaf layer
end

%-----------------------------------------------------------------------
% --- Longwave radiation transfer through canopy (analytical method)
%-----------------------------------------------------------------------

% Diffuse (Kd) and direct beam (Kb) extinction coefficients

Kd = 0.78;
Kb = 0.5;

% Longwave flux from ground and leaf

Lgrnd = emgrnd * sigma * tgrnd^4;
Lleaf = emleaf * sigma * tveg^4;

% Canopy integration: compare analytical solution with numerical integration

Lc = emleaf * (irsky + Lgrnd) * (1 - exp(-Kd*LAI)) - 2 * Lleaf * (1 - exp(-Kd*LAI));

f1 = @(x) (emleaf*Lgrnd - Lleaf) * Kd * exp(-Kd * (LAI-x)) + (emleaf*irsky - Lleaf) * Kd * exp(-Kd * x);
Lc_numerical = integral(f1, 0, LAI);

fprintf('Analytical model \n')
fprintf('Lc = %15.5f\n',Lc)
fprintf('Lc = %15.5f\n',Lc_numerical)

irveg = Lc;

% Sunlit canopy: compare analytical solution with numerical integration

Lcsun = (emleaf * irsky - Lleaf) * Kd / (Kd + Kb) * (1 - exp(-(Kd+Kb)*LAI)) + ...
   (emleaf * Lgrnd - Lleaf) * Kd / (Kd - Kb) * (exp(-Kb*LAI) - exp(-Kd*LAI));

f1sun = @(x) f1(x) .* exp(-Kb * x);
Lcsun_numerical = integral(f1sun, 0, LAI);

fprintf('Lcsun = %15.5f\n',Lcsun)
fprintf('Lcsun = %15.5f\n',Lcsun_numerical)

% Shaded canopy: compare analytical solution with numerical integration

Lcsha = Lc - Lcsun;

f1sha = @(x) f1(x) .* (1 - exp(-Kb * x));
Lcsha_numerical = integral(f1sha, 0, LAI);

fprintf('Lcsha = %15.5f\n',Lcsha)
fprintf('Lcsha = %15.5f\n',Lcsha_numerical)

% Absorbed longwave radiation for ground (soil)

Ld = irsky * (1 - emleaf * (1 - exp(-Kd * LAI))) + emleaf * sigma * tveg^4 * (1 - exp(-Kd * LAI));
irsoi = Ld - Lgrnd;

% Canopy emitted longwave radiation

Lu = Lgrnd * (1 - emleaf * (1 - exp(-Kd * LAI))) + emleaf * sigma * tveg^4 * (1 - exp(-Kd * LAI));
irup = Lu;

% Conservation check: absorbed = incoming - outgoing

sumabs = irsky - irup;
err = sumabs - (irveg + irsoi);
if (abs(err) > 1e-03)
   fprintf('err = %15.5f\n',err)
   fprintf('sumabs = %15.5f\n',sumabs)
   fprintf('irveg = %15.5f\n',irveg)
   fprintf('irsoi = %15.5f\n',irsoi)
   error ('Analytical solution: Longwave conservation error')
end

fprintf(' \n')
fprintf('irup = %15.5f\n',irup)
fprintf('irveg = %15.5f\n',irveg)
fprintf('irsoi = %15.5f\n',irsoi)
fprintf(' \n')

%-----------------------------------------------------------------------
% --- Longwave radiation transfer through canopy using Norman (1979)
%-----------------------------------------------------------------------

fprintf('Norman model \n')

% --- Leaf scattering coefficient

omega = 1 - emleaf;

% --- Intercepted radiation is reflected

rho = omega;   % Leaf reflectance
tau = 0;       % Leaf transmittance

% --- Intercepted radiation is both reflected and transmitted

% rho = omega * 0.5;
% tau = omega * 0.5;

% --- Emitted longwave radiation from leaves (W/m2)

for iv = nbot:ntop
   ir_source(iv) = emleaf * sigma * tleaf(iv)^4 * (1 - td(iv));
end

% --- Set up and solve tridiagonal system of equations for upward and downward fluxes

% There are two equations for each canopy layer and the soil. The first
% equation is the upward flux and the second equation is the downward flux.

m = 0;

% Soil: upward flux

iv = nsoi;
m = m + 1;
atri(m) = 0;
btri(m) = 1;
ctri(m) = -(1 - emgrnd);
dtri(m) = emgrnd * sigma * tgrnd^4;

% Soil: downward flux

refld = (1 - td(iv+1)) * rho;
trand = (1 - td(iv+1)) * tau + td(iv+1);
aiv = refld - trand * trand / refld;
biv = trand / refld;

m = m + 1;
atri(m) = -aiv;
btri(m) = 1;
ctri(m) = -biv;
dtri(m) = (1 - biv) * ir_source(iv+1);

% Leaf layers, excluding top layer

for iv = nbot:ntop-1

   % Upward flux

   refld = (1 - td(iv)) * rho;
   trand = (1 - td(iv)) * tau + td(iv);
   fiv = refld - trand * trand / refld;
   eiv = trand / refld;

   m = m + 1;
   atri(m) = -eiv;
   btri(m) = 1;
   ctri(m) = -fiv;
   dtri(m) = (1 - eiv) * ir_source(iv);

   % Downward flux

   refld = (1 - td(iv+1)) * rho;
   trand = (1 - td(iv+1)) * tau + td(iv+1);
   aiv = refld - trand * trand / refld;
   biv = trand / refld;

   m = m + 1;
   atri(m) = -aiv;
   btri(m) = 1;
   ctri(m) = -biv;
   dtri(m) = (1 - biv) * ir_source(iv+1);

end

% Top canopy layer: upward flux

iv = ntop;
refld = (1 - td(iv)) * rho;
trand = (1 - td(iv)) * tau + td(iv);
fiv = refld - trand * trand / refld;
eiv = trand / refld;

m = m + 1;
atri(m) = -eiv;
btri(m) = 1;
ctri(m) = -fiv;
dtri(m) = (1 - eiv) * ir_source(iv);

% Top canopy layer: downward flux

m = m + 1;
atri(m) = 0;
btri(m) = 1;
ctri(m) = 0;
dtri(m) = irsky;

% Solve tridiagonal equations for upward and downward fluxes

[utri] = tridiagonal_solver (atri, btri, ctri, dtri, m);

% Now copy the solution (utri) to the upward (irup) and downward (irdn)
% fluxes for each layer
% irup - Upward longwave flux above layer
% irdn - Downward longwave flux onto layer

m = 0;

% Soil fluxes

iv = nsoi;
m = m + 1;
irup(iv) = utri(m);
m = m + 1;
irdn(iv) = utri(m);

% Leaf layer fluxes

for iv = nbot:ntop
   m = m + 1;
   irup(iv) = utri(m);
    m = m + 1;
   irdn(iv) = utri(m);
end

% --- Error check: compare tridiagonal solution with actual equations

iv = ntop;
irdn_eq(iv) = irsky;
for iv = ntop-1: -1: nsoi
   irdn_eq(iv) = irdn(iv+1) * (td(iv+1)+(1-td(iv+1))*tau) + irup(iv) * ((1-td(iv+1))*rho) ...
               + emleaf*sigma*tleaf(iv+1)^4*(1-td(iv+1));
end

iv = nsoi;
irup_eq(iv) = (1-emgrnd) * irdn(iv) + emgrnd*sigma*tgrnd^4;
for iv = nsoi:ntop-1
   irup_eq(iv+1) = irup(iv) * (td(iv+1)+(1-td(iv+1))*tau) + irdn(iv+1) * ((1-td(iv+1))*rho) ...
                 + emleaf*sigma*tleaf(iv+1)^4*(1-td(iv+1));
end

for iv = nsoi:ntop
   err = irdn_eq(iv) - irdn(iv);
   if (abs(err) > 1e-10)
      fprintf('err = %15.5f\n',err)
      fprintf('tridiag = %15.5f\n',irdn(iv))
      fprintf('eq = %15.5f\n',irdn_eq(iv))
      error ('Norman radiation: downward error')
   end
   err = irup_eq(iv) - irup(iv);
   if (abs(err) > 1e-10)
      fprintf('err = %15.5f\n',err)
      fprintf('tridiag = %15.5f\n',irup(iv))
      fprintf('eq = %15.5f\n',irup_eq(iv))
      error ('Norman radiation: upward error')
   end
end

% --- Compute fluxes

% Absorbed longwave radiation for ground (soil)

iv = nsoi;
irabs(iv) = irdn(iv) - irup(iv);
irsoi = irabs(iv);

% Absorbed longwave radiation for leaf layers

for iv = nbot:ntop
   irabs(iv) = emleaf * (irdn(iv) + irup(iv-1)) * (1 - td(iv)) - 2 * ir_source(iv);
end

% Sum longwave radiation absorbed by vegetation

irveg = 0;
for iv = nbot:ntop
   irveg = irveg + irabs(iv);
end

% Canopy emitted longwave radiation

irup = irup(ntop);

% --- Conservation check: absorbed = incoming - outgoing

sumabs = irsky - irup;
err = sumabs - (irveg + irsoi);
if (abs(err) > 1e-03)
   fprintf('err = %15.5f\n',err)
   fprintf('sumabs = %15.5f\n',sumabs)
   fprintf('irveg = %15.5f\n',irveg)
   fprintf('irsoi = %15.5f\n',irsoi)
   error ('NormanRadiation: Longwave conservation error')
end

fprintf('irup = %15.5f\n',irup)
fprintf('irveg = %15.5f\n',irveg)
fprintf('irsoi = %15.5f\n',irsoi)

function [u] = tridiagonal_solver (a, b, c, d, n)

% Solve for U given the set of equations R * U = D, where U is a vector
% of length N, D is a vector of length N, and R is an N x N tridiagonal
% matrix defined by the vectors A, B, C each of length N. A(1) and
% C(N) are undefined and are not referenced.
%
%     |B(1) C(1) ...  ...  ...                     |
%     |A(2) B(2) C(2) ...  ...                     |
% R = |     A(3) B(3) C(3) ...                     |
%     |                    ... A(N-1) B(N-1) C(N-1)|
%     |                    ... ...    A(N)   B(N)  |
%
% The system of equations is written as:
%
%    A_i * U_i-1 + B_i * U_i + C_i * U_i+1 = D_i
%
% for i = 1 to N. The solution is found by rewriting the
% equations so that:
%
%    U_i = F_i - E_i * U_i+1

% --- Forward sweep (1 -> N) to get E and F

e(1) = c(1) / b(1);

for i = 2: 1: n-1
   e(i) = c(i) / (b(i) - a(i) * e(i-1));
end

f(1) = d(1) / b(1);

for i = 2: 1: n
   f(i) = (d(i) - a(i) * f(i-1)) / (b(i) - a(i) * e(i-1));
end

% --- Backward substitution (N -> 1) to solve for U

u(n) = f(n);

for i = n-1: -1: 1
   u(i) = f(i) - e(i) * u(i+1);
end

end
