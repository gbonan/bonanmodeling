function [psi] = matric_potential (type, params, theta)

% --- Calculate psi for a given theta

switch type
   case 'van_Genuchten'

   theta_res = params(1);    % Residual water content
   theta_sat = params(2);    % Volumetric water content at saturation
   alpha = params(3);        % Inverse of the air entry potential
   n = params(4);            % Pore-size distribution index
   m = params(5);            % Exponent

   Se = (theta - theta_res) / (theta_sat - theta_res);
   psi = -((Se^(-1/m) - 1)^(1/n)) / alpha;

   case 'Campbell'

   theta_sat = params(1);    % Volumetric water content at saturation
   psi_sat = params(2);      % Matric potential at saturation
   b = params(3);            % Exponent

   psi = psi_sat * (theta / theta_sat)^-b;

end
