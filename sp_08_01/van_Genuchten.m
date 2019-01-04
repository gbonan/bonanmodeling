function [theta, K, cap] = van_Genuchten (params, psi)

% ----------------------------------
% van Genuchten (1980) relationships
% ----------------------------------

% --- Soil parameters

theta_res = params(1);   % Residual water content
theta_sat = params(2);   % Volumetric water content at saturation
alpha = params(3);       % Inverse of the air entry potential
n = params(4);           % Pore-size distribution index
m = params(5);           % Exponent
Ksat = params(6);        % Hydraulic conductivity at saturation
ityp = params(7);        % Soil texture flag

% --- Effective saturation (Se) for specified matric potential (psi)

if (psi <= 0)
   Se = (1 + (alpha * abs(psi))^n)^-m;
else
   Se = 1;
end

% --- Volumetric soil moisture (theta) for specified matric potential (psi)

theta = theta_res + (theta_sat - theta_res) * Se;

% --- Hydraulic conductivity (K) for specified matric potential (psi)

if (Se <= 1)
   K = Ksat * sqrt(Se) * (1 - (1 - Se^(1/m))^m)^2;

   % Special case for Haverkamp et al. (1977) sand (ityp = 1) and Yolo light clay (ityp = 2)

   if (ityp == 1)
      K = Ksat * 1.175e6 / (1.175e6 + abs(psi)^4.74);
   end
   if (ityp == 2)
      K = Ksat * 124.6/ (124.6 + abs(psi)^1.77);
   end

else

   K = Ksat;

end

% --- Specific moisture capacity (cap) for specified matric potential (psi)

if (psi <= 0)
   num = alpha * m * n * (theta_sat - theta_res) * (alpha * abs(psi))^(n-1);
   den =  (1 + (alpha * abs(psi))^n)^(m+1);
   cap = num / den;
else
   cap = 0;
end
