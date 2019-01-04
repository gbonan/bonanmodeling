function [val] = latvap (tc, mmh2o)

% Latent heat of vaporization (J/mol) at temperature tc (degC)

val = 2501.6 - 2.3773 * tc; % Latent heat of vaporization (J/g)
val = val * 1000;           % Convert from J/g to J/kg
val = val * mmh2o;          % Convert from J/kg to J/mol
