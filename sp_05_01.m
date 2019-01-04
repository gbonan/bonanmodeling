% Supplemental program 5.1

% ----------------------------------------
% Calculate and graph thermal conductivity
% ----------------------------------------

cwat = 4188;               % Specific heat of water (J/kg/K)
cice = 2117.27;            % Specific heat ice (J/kg/K)

rho_wat = 1000;            % Density of water (kg/m3)
rho_ice = 917;             % Density of ice (kg/m3)

cvwat = cwat * rho_wat;    % Heat capacity of water (J/m3/K)
cvice = cice * rho_ice;    % Heat capacity of ice (J/m3/K)
cvsol = 1.926e06;          % Heat capacity of soil solids (J/m3/K)

tkwat = 0.57;              % Thermal conductivity of water (W/m/K)
tkice = 2.29;              % Thermal conductivity of ice (W/m/K)
tk_quartz = 7.7;           % Thermal conductivity of quartz (W/m/K)

% Soil texture classes (Cosby et al. 1984. Water Resources Research 20:682-690)

%  1: sand
%  2: loamy sand
%  3: sandy loam
%  4: silty loam
%  5: loam
%  6: sandy clay loam
%  7  silty clay loam
%  8: clay loam
%  9: sandy clay
% 10: silty clay
% 11: clay

silt = [ 5, 12, 32, 70, 39, 15, 56, 34,  6, 47, 20];       % Percent silt
sand = [92, 82, 58, 17, 43, 58, 10, 32, 52,  6, 22];       % Percent sand
clay = [ 3,  6, 10, 13, 18, 27, 34, 34, 42, 47, 58];       % Percent clay

% Volumetric soil water content at saturation (porosity)
% (Clapp and Hornberger. 1978. Water Resources Research 14:601-604)

watsat = [0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.492, 0.482];

% Define 5 soil types to process

soiltyp = [1, 3, 5, 8, 11];

% Set relative soil water content (s) from 0 to 1

inc = 0.05;                             % increment
n = (1 - 0) / inc + 1;                  % number of values
s = linspace(0,1,n);                    % n evenly spaced values between 0 and 1 (inclusive)

% Loop through each soil type

for i = 1:length(soiltyp)

   % Soil texture to process

   k = soiltyp(i);

   % Thermal conductivity and heat capacity for each soil moisture

   for j = 1:length(s)

      % Volumetric water content

      h2osoi = s(j) * watsat(k);

      % Dry thermal conductivity (W/m/K) from bulk density (kg/m3)

      bd = 2700 * (1 - watsat(k));
      tkdry = (0.135 * bd + 64.7) / (2700 - 0.947 * bd);

      % Soil solids thermal conducitivty (W/m/K) from quartz fraction
      % tko = thermal conductivity of other minerals (W/m/K)

      quartz = sand(k) / 100;
      if (quartz > 0.2)
         tko = 2.0;
      else
         tko = 3.0;
      end
      tksol = (tk_quartz^quartz) * (tko^(1-quartz));

      % Unfrozen and frozen saturated thermal conductivity (W/m/K)

      tksat_u = (tksol^(1-watsat(k))) * (tkwat^watsat(k));
      tksat_f = (tksol^(1-watsat(k))) * (tkice^watsat(k));

      % Unfrozen and frozen Kersten number

      if (sand(k) < 50)
         ke_u = log10(max(s(j),0.1)) + 1;
      else
         ke_u = 0.7 * log10(max(s(j),0.05)) + 1;
      end
      ke_f = s(j);

      % Unfrozen and frozen thermal conductivity (W/m/K)

      tku = (tksat_u - tkdry) * ke_u + tkdry;
      tkf = (tksat_f - tkdry) * ke_f + tkdry;

      % Unfrozen and frozen heat capacity (J/m3/K)

      cvu = (1 - watsat(k)) * cvsol + cvwat * h2osoi;
      cvf = (1 - watsat(k)) * cvsol + cvice * h2osoi;

      % Save values for each texture type

      if (i == 1)
         tk1(j) = tku;
         cv1(j) = cvu * 1e-06;
      elseif (i == 2)
         tk2(j) = tku;
         cv2(j) = cvu * 1e-06;
      elseif (i == 3)
         tk3(j) = tku;
         cv3(j) = cvu * 1e-06;
      elseif (i == 4)
         tk4(j) = tku;
         cv4(j) = cvu * 1e-06;
      elseif (i == 5)
         tk5(j) = tku;
         cv5(j) = cvu * 1e-06;
      end

   end      % end soil water loop j
end         % end soil type loop i

% Make graph

plot(s,tk1,'r-',s,tk2,'g-',s,tk3,'b-',s,tk4,'m-',s,tk5,'c-')
title('Thermal conductivity')
xlabel('Relative soil moisture')
ylabel('Thermal conductivity (W m^{-1} K^{-1})')
legend('sand','sandy loam','loam','clay loam','clay','Location','best')

% Write formated output to text file: n rows x 6 columns
% column 1 = relative soil water (s)
% column 2 = thermal conductivity for soil 1 (W/m/K)
% column 3 = thermal conductivity for soil 2 (W/m/K)
% column 4 = thermal conductivity for soil 3 (W/m/K)
% column 5 = thermal conductivity for soil 4 (W/m/K)
% column 6 = thermal conductivity for soil 5 (W/m/K)

A = [s; tk1; tk2; tk3; tk4; tk5];
fileID = fopen('tk.txt','w');
fprintf(fileID,'%12s %12s %12s %12s %12s %12s\n','s','tex1','tex2','tex3','tex4','tex5');
fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n', A);
fclose(fileID);

% Write formated output to text file: n rows x 6 columns
% column 1 = relative soil water (s)
% column 2 = heat capacity for soil 1 (MJ/m3/K)
% column 3 = heat capacity for soil 2 (MJ/m3/K)
% column 4 = heat capacity for soil 3 (MJ/m3/K)
% column 5 = heat capacity for soil 4 (MJ/m3/K)
% column 6 = heat capacity for soil 5 (MJ/m3/K)

B = [s; cv1; cv2; cv3; cv4; cv5];
fileID = fopen('cv.txt','w');
fprintf(fileID,'%12s %12s %12s %12s %12s %12s\n','s','tex1','tex2','tex3','tex4','tex5');
fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n', B);
fclose(fileID);
