% Supplemental program 6.1

% -------------------------------------------------------------------------
% Calculate friction velocity and sensible heat flux given wind speed and
% temperature at two heights using Monin-Obukhov similarity theory or
% roughness sublayer theory from Physick and Garratt (1995)
% -------------------------------------------------------------------------

% --- Physical constants

rgas = 8.31446;           % Universal gas constant (J/K/mol)
var.k = 0.4;              % von Karman constant
var.g = 9.80665;          % Gravitational acceleration (m/s2)
cpair = 29.2;             % Specific heat of air at constant pressure (J/mol/K)

% --- Input variables

var.d = 19.0;             % Displacement height (m)

var.z1 = 21.0;            % Height (m)
var.u1 = 1.0;             % Wind speed at height z1 (m/s)
var.t1 = 29.0 + 273.15;   % Temperature at height z1 (K)

var.z2 = 29.0;            % Height (m)
var.u2 = 2.1;             % Wind speed at height z2 (m/s)
var.t2 = 28.1 + 273.15;   % Temperature at height z2 (K)

var.zstar = 49.0;         % Height of roughness sublayer (m)

  abl = 'MOST';           % Use Monin-Obukhov similarity theory
% abl = 'RSL';            % Use roughness sublayer theory

% --- Molar density (mol/m3)

rhomol = 101325 / (rgas * var.t2);

switch abl

   % -----------------------------------
   % Monin-Obukhov similarity theory
   % -----------------------------------

   case 'MOST'

   % Use bisection to solve for L as specified by the function "most"
   % and then calculate fluxes for that value of L

   L1 = 100;                                    % Initial guess for Obukhov length L (m)
   L2 = -100;                                   % Initial guess for Obukhov length L (m)
   func_name = 'most';                          % The function name is "most", in the file most.m
   [L] = bisect (func_name, L1, L2, 0.01, var); % Solve for L (m)

   % Evaluate psi for momentum and scalars at heights z2 and z1

   [psi_m_z2] = psi_m_monin_obukhov((var.z2-var.d)/L);
   [psi_m_z1] = psi_m_monin_obukhov((var.z1-var.d)/L);
   [psi_c_z2] = psi_c_monin_obukhov((var.z2-var.d)/L);
   [psi_c_z1] = psi_c_monin_obukhov((var.z1-var.d)/L);

   % Calculate u* and T* and the sensible heat flux

   ustar = (var.u2 - var.u1) * var.k / (log((var.z2-var.d)/(var.z1-var.d)) - (psi_m_z2 - psi_m_z1));
   tstar = (var.t2 - var.t1) * var.k / (log((var.z2-var.d)/(var.z1-var.d)) - (psi_c_z2 - psi_c_z1));
   H = -rhomol * cpair * tstar * ustar;

   % Calculate aerodynamic conductances

   gam = rhomol * var.k * ustar / (log((var.z2-var.d)/(var.z1-var.d)) - (psi_m_z2 - psi_m_z1));
   gac = rhomol * var.k * ustar / (log((var.z2-var.d)/(var.z1-var.d)) - (psi_c_z2 - psi_c_z1));

   fprintf('Monin-Obuhkov similarity theory\n')
   fprintf('L = %15.3f\n',L)
   fprintf('u* = %15.3f\n',ustar)
   fprintf('T* = %15.3f\n',tstar)
   fprintf('H = %15.3f\n',H)
   fprintf('gam = %15.3f\n',gam)
   fprintf('gac = %15.3f\n',gac)

   % -----------------------------------
   % Roughness sublayer theory
   % -----------------------------------

   case 'RSL'

   % Use bisection to solve for L as specified by the function "rsl"
   % and then calculate fluxes for that value of L

   L1 = 100;                                    % Initial guess for Obukhov length L (m)
   L2 = -100;                                   % Initial guess for Obukhov length L (m)
   func_name = 'rsl';                           % The function name is "rsl", in the file rsl.m
   [L] = bisect (func_name, L1, L2, 0.01, var); % Solve for L (m)

   % Evaluate psi for momentum and scalars at heights z2 and z1

   [psi_m_z2] = psi_m_monin_obukhov((var.z2-var.d)/L);
   [psi_m_z1] = psi_m_monin_obukhov((var.z1-var.d)/L);
   [psi_c_z2] = psi_c_monin_obukhov((var.z2-var.d)/L);
   [psi_c_z1] = psi_c_monin_obukhov((var.z1-var.d)/L);

   % Evaluate the roughness sublayer-modified psi (between z1 and z2)

   f1_psi_m_rsl = @(z) (1-16*(z-var.d)/L).^(-0.25) .* (1-exp(-0.7*(1-(z-var.d)/(var.zstar-var.d)))) ./ (z-var.d);
   f1_psi_c_rsl = @(z) (1-16*(z-var.d)/L).^(-0.50) .* (1-exp(-0.7*(1-(z-var.d)/(var.zstar-var.d)))) ./ (z-var.d);

   f2_psi_m_rsl = @(z) (1+5*(z-var.d)/L) .* (1-exp(-0.7*(1-(z-var.d)/(var.zstar-var.d)))) ./ (z-var.d);
   f2_psi_c_rsl = @(z) (1+5*(z-var.d)/L) .* (1-exp(-0.7*(1-(z-var.d)/(var.zstar-var.d)))) ./ (z-var.d);

   if (L < 0)
      psi_m_rsl = integral (f1_psi_m_rsl, var.z1, var.z2);
      psi_c_rsl = integral (f1_psi_c_rsl, var.z1, var.z2);
   else
      psi_m_rsl = integral (f2_psi_m_rsl, var.z1, var.z2);
      psi_c_rsl = integral (f2_psi_c_rsl, var.z1, var.z2);
   end

   % Calculate u* and T* and the sensible heat flux

   ustar = (var.u2 - var.u1) * var.k / (log((var.z2-var.d)/(var.z1-var.d)) - (psi_m_z2 - psi_m_z1) - psi_m_rsl);
   tstar = (var.t2 - var.t1) * var.k / (log((var.z2-var.d)/(var.z1-var.d)) - (psi_c_z2 - psi_c_z1) - psi_c_rsl);
   H = -rhomol * cpair * tstar * ustar;

   % Calculate aerodynamic conductances

   gam = rhomol * var.k * ustar / (log((var.z2-var.d)/(var.z1-var.d)) - (psi_m_z2 - psi_m_z1) - psi_m_rsl);
   gac = rhomol * var.k * ustar / (log((var.z2-var.d)/(var.z1-var.d)) - (psi_c_z2 - psi_c_z1) - psi_c_rsl);

   fprintf('Roughness sublayer theory\n')
   fprintf('L = %15.3f\n',L)
   fprintf('u* = %15.3f\n',ustar)
   fprintf('T* = %15.3f\n',tstar)
   fprintf('H = %15.3f\n',H)
   fprintf('gam = %15.3f\n',gam)
   fprintf('gac = %15.3f\n',gac)

end
