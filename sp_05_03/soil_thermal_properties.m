function [soilvar] = soil_thermal_properties (physcon, soilvar)

% Calculate soil thermal conductivity and heat capacity

% ------------------------------------------------------
% Input
%   physcon.hfus             ! Heat of fusion for water at 0 C (J/kg)
%   physcon.tfrz             ! Freezing point of water (K)
%   physcon.rhowat           ! Density of water (kg/m3)
%   physcon.rhoice           ! Density of ice (kg/m3)
%   soilvar.method           ! Use excess heat or apparent heat capacity for phase change
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

% Temperature range for freezing and thawing (K)

tinc = 0.5;

% Unfrozen and frozen thermal conductivity (W/m/K)

tku = 1.860;
tkf = 2.324;

% Unfrozen and frozen heat capacity (J/m3/K)

cvu = 2.862e06;
cvf = 1.966e06;

for i = 1:soilvar.nsoi

   % --- Volumetric soil water and ice

   watliq = soilvar.h2osoi_liq(i) / (physcon.rhowat * soilvar.dz(i));
   watice = soilvar.h2osoi_ice(i) / (physcon.rhoice * soilvar.dz(i));

   % Heat of fusion (J/m3) - This is equivalent to ql = hfus * (h2osoi_liq + h2osoi_ice) / dz

   ql = physcon.hfus * (physcon.rhowat * watliq + physcon.rhoice * watice);

   % Heat capacity and thermal conductivity

   if (soilvar.tsoi(i) > physcon.tfrz+tinc)
      soilvar.cv(i) = cvu;
      soilvar.tk(i) = tku;
   end

   if (soilvar.tsoi(i) >= physcon.tfrz-tinc & soilvar.tsoi(i) <= physcon.tfrz+tinc)
      switch soilvar.method
         case 'apparent-heat-capacity'
         soilvar.cv(i) = (cvf + cvu) / 2 + ql / (2 * tinc);
         case 'excess heat'
         soilvar.cv(i) = (cvf + cvu) / 2;
      end
      soilvar.tk(i) = tkf + (tku - tkf) * (soilvar.tsoi(i) - physcon.tfrz + tinc) / (2 * tinc);
   end

   if (soilvar.tsoi(i) < physcon.tfrz-tinc)
      soilvar.cv(i) = cvf;
      soilvar.tk(i) = tkf;
   end

end
