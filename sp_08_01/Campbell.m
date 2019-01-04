function [theta, K, cap] = Campbell (params, psi)

% -----------------------------
% Campbell (1974) relationships
% -----------------------------

% --- Soil parameters

theta_sat = params(1);    % Volumetric water content at saturation
psi_sat = params(2);      % Matric potential at saturation
b = params(3);            % Exponent
Ksat = params(4);         % Hydraulic conductivity at saturation

% --- Volumetric soil moisture (theta) for specified matric potential (psi)

if (psi <= psi_sat)
   theta = theta_sat * (psi / psi_sat)^(-1/b);
else
   theta = theta_sat;
end

% --- Hydraulic conductivity (K) for specified matric potential (psi)

if (psi <= psi_sat)
   K = Ksat * (theta / theta_sat)^(2*b+3);
else
   K = Ksat;
end

% --- Specific moisture capacity (cap) for specified matric potential (psi)

if (psi <= psi_sat)
   cap = -theta_sat / (b * psi_sat) * (psi / psi_sat)^(-1/b-1);
else
   cap = 0;
end
