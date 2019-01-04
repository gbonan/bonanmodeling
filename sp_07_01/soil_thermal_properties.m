function [soilvar] = soil_thermal_properties (physcon, soilvar)

% Calculate soil thermal conductivity and heat capacity

% ------------------------------------------------------
% Input
%   physcon.hfus             ! Heat of fusion for water at 0 C (J/kg)
%   physcon.tfrz             ! Freezing point of water (K)
%   physcon.tkwat            ! Thermal conductivity of water (W/m/K)
%   physcon.tkice            ! Thermal conductivity of ice (W/m/K)
%   physcon.cvwat            ! Heat capacity of water (J/m3/K)
%   physcon.cvice            ! Heat capacity of ice (J/m3/K)
%   physcon.rhowat           ! Density of water (kg/m3)
%   physcon.rhoice           ! Density of ice (kg/m3)
%   soilvar.method           ! Use excess heat or apparent heat capacity for phase change
%   soilvar.soil_texture     ! Soil texture class
%   soilvar.sand             ! Percent sand
%   soilvar.watsat           ! Volumetric soil water content at saturation (porosity)
%   soilvar.nsoi             ! Number of soil layers
%   soilvar.dz               ! Soil layer thickness (m)
%   soilvar.tsoi             ! Soil temperature (K)
%   soilvar.h2osoi_liq       ! Unfrozen water, liquid (kg H2O/m2)
%   soilvar.h2osoi_ice       ! Frozen water, ice (kg H2O/m2)
%
% Input/output
%   soilvar.tk               ! Thermal conducitivty (W/m/K)
%   soilvar.cv               ! Volumetric heat capacity (J/m3/K)
% ------------------------------------------------------

for i = 1:soilvar.nsoi

   % --- Soil texture to process

   k = soilvar.soil_texture;

   % --- Volumetric soil water and ice

   watliq = soilvar.h2osoi_liq(i) / (physcon.rhowat * soilvar.dz(i));
   watice = soilvar.h2osoi_ice(i) / (physcon.rhoice * soilvar.dz(i));

   % Fraction of total volume that is liquid water

   fliq = watliq / (watliq + watice);

   % Soil water relative to saturation

   s = min((watliq + watice)/soilvar.watsat(k), 1);

   % --- Dry thermal conductivity (W/m/K) from bulk density (kg/m3)

   bd = 2700 * (1 - soilvar.watsat(k));
   tkdry = (0.135 * bd + 64.7) / (2700 - 0.947 * bd);

   % --- Soil solids thermal conducitivty (W/m/K)

   % Thermal conductivity of quartz (W/m/K)

   tk_quartz = 7.7;

   % Quartz fraction

   quartz = soilvar.sand(k) / 100;

   % Thermal conductivity of other minerals (W/m/K)

   if (quartz > 0.2)
      tko = 2;
   else
      tko = 3;
   end

   % Thermal conductivity of soil solids (W/m/K)

   tksol = tk_quartz^quartz * tko^(1-quartz);

   % --- Saturated thermal conductivity (W/m/K) and unfrozen and frozen values

   tksat = tksol^(1-soilvar.watsat(k)) * physcon.tkwat^(fliq*soilvar.watsat(k)) * ...
      physcon.tkice^(soilvar.watsat(k)-fliq*soilvar.watsat(k));
   tksat_u = tksol^(1-soilvar.watsat(k)) * physcon.tkwat^soilvar.watsat(k);
   tksat_f = tksol^(1-soilvar.watsat(k)) * physcon.tkice^soilvar.watsat(k);

   % --- Kersten number and unfrozen and frozen values

   if (soilvar.sand(k) < 50)
      ke_u = log10(max(s,0.1)) + 1;
   else
      ke_u = 0.7 * log10(max(s,0.05)) + 1;
   end
   ke_f = s;

   if (soilvar.tsoi(i) >= physcon.tfrz)
      ke = ke_u;
   else
      ke = ke_f;
   end

   % --- Thermal conductivity (W/m/K) and unfrozen and frozen values

   soilvar.tk(i) = (tksat - tkdry) * ke + tkdry;
   tku = (tksat_u - tkdry) * ke_u + tkdry;
   tkf = (tksat_f - tkdry) * ke_f + tkdry;

   % --- Heat capacity of soil solids (J/m3/K)

   cvsol = 1.926e06;

   % --- Heat capacity (J/m3/K) and unfrozen and frozen values

   soilvar.cv(i) = (1 - soilvar.watsat(k)) * cvsol + physcon.cvwat * watliq + physcon.cvice * watice;
   cvu = (1 - soilvar.watsat(k)) * cvsol + physcon.cvwat * (watliq + watice);
   cvf = (1 - soilvar.watsat(k)) * cvsol + physcon.cvice * (watliq + watice);

   % --- Adjust heat capacity and thermal conductivity if using apparent heat capacity

   switch soilvar.method
      case 'apparent-heat-capacity'

      % Temperature range for freezing and thawing (K)

      tinc = 0.5;

      % Heat of fusion (J/m3) - This is equivalent to ql = hfus * (h2osoi_liq + h2osoi_ice) / dz

      ql = physcon.hfus * (physcon.rhowat * watliq + physcon.rhoice * watice);

      % Heat capacity and thermal conductivity

      if (soilvar.tsoi(i) > physcon.tfrz+tinc)
         soilvar.cv(i) = cvu;
         soilvar.tk(i) = tku;
      end

      if (soilvar.tsoi(i) >= physcon.tfrz-tinc & soilvar.tsoi(i) <= physcon.tfrz+tinc)
         soilvar.cv(i) = (cvf + cvu) / 2 + ql / (2 * tinc);
         soilvar.tk(i) = tkf + (tku - tkf) * (soilvar.tsoi(i) - physcon.tfrz + tinc) / (2 * tinc);
      end

      if (soilvar.tsoi(i) < physcon.tfrz-tinc)
         soilvar.cv(i) = cvf;
         soilvar.tk(i) = tkf;
      end
   end

end
