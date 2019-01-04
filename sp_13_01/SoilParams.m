function [soil] = SoilParams (soil)

% Set soil hydraulic parameters, soil depth, and rooting fraction

% -------------------------------------------------------------------------
% Input
%   soil.texture     ! Soil texture class
% Output
%   soil.nlevsoi     ! Number of soil layers
%   soil.dz          ! Soil layer thickness (m)
%   soil.rootfr      ! Fraction of roots in each soil layer (-)
%   soil.watsat      ! Soil layer volumetric water content at saturation (porosity)
%   soil.psisat      ! Soil layer matric potential at saturation (mm)
%   soil.hksat       ! Soil layer hydraulic conductivity at saturation (mm H2O/s)
%   soil.bsw         ! Soil layer Clapp and Hornberger "b" parameter
% -------------------------------------------------------------------------

% --- Soil texture classes
%  1: sand
%  2: loamy sand
%  3: sandy loam
%  4: silt loam
%  5: loam
%  6: sandy clay loam
%  7  silty clay loam
%  8: clay loam
%  9: sandy clay
% 10: silty clay
% 11: clay

% Volumetric water content at saturation (porosity)

theta_sat = [0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.492, 0.482];

% Matric potential at saturation (mm)

psi_sat = [-121, -90, -218, -786, -478, -299, -356, -630, -153, -490, -405];

% Clapp and Hornberger "b" parameter

b = [4.05, 4.38, 4.90, 5.30, 5.39, 7.12, 7.75, 8.52, 10.40, 10.40, 11.40];

% Hydraulic conductivity at saturation (cm/min -> mm/s)

k_sat = [1.056, 0.938, 0.208, 0.0432, 0.0417, 0.0378, 0.0102, 0.0147, 0.0130, 0.0062, 0.0077];
k_sat = k_sat * 10 / 60;

% --- Soil layer thickness (total depth is 2.5 m)

soil.nlevsoi = 11;
soil.dz = [0.05, 0.05, 0.10, 0.10, 0.20, 0.20, 0.20, 0.30, 0.40, 0.40, 0.50];

% --- Root profile

  beta_param = 0.90;  % root profile parameter: shallow profile
% beta_param = 0.97;  % root profile parameter: deep profile

for j = 1:soil.nlevsoi
   if (j == 1)
      z2 = soil.dz(j) * 100;
      soil.rootfr(j) = 1 - beta_param^z2;
   else
      z1 = z2;
      z2 = z1 + soil.dz(j) * 100;
      soil.rootfr(j) = beta_param^z1 - beta_param^z2;
   end
end

% --- Soil texture to process

k = soil.texture;

% --- Soil layer values

for j = 1:soil.nlevsoi
   soil.watsat(j) = theta_sat(k);
   soil.psisat(j) = psi_sat(k);
   soil.hksat(j) = k_sat(k);
   soil.bsw(j) = b(k);
end
