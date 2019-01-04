% Supplemental program 12.1

% -------------------------------------------------------------------------
% Calculate leaf gas exchange for C3 and C4 plants in relation to PAR, CO2,
% temperature, and vapor pressure deficit. For this example, leaf temperature
% is specified (not calculated).
% -------------------------------------------------------------------------

% --- Waveband indices for visible and near-infrared

params.vis = 1; params.nir = 2;

% --- Physical constants

physcon.grav = 9.80665;               % Gravitational acceleration (m/s2)
physcon.tfrz = 273.15;                % Freezing point of water (K)
physcon.mmh2o = 18.02 / 1000;         % Molecular mass of water (kg/mol)
physcon.rgas = 8.31446;               % Universal gas constant (J/K/mol)
physcon.visc0 = 13.3e-06;             % Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
physcon.Dh0 = 18.9e-06;               % Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
physcon.Dv0 = 21.8e-06;               % Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
physcon.Dc0 = 13.8e-06;               % Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)

% --- Set leaf physiology variables

% Stomatal conductance: 0 = Medlyn model. 1 = Ball-Berry model. 2 = WUE optimization

  leaf.gstyp = 0;
% leaf.gstyp = 1;
% leaf.gstyp = 2;

% Photosynthetic pathway: 1 = C3. 0 = C4

  leaf.c3psn = 1;
% leaf.c3psn = 0;

% Photosynthesis co-limitation: 0 = no. 1 = yes

  leaf.colim = 1;
% leaf.colim = 0;

% Leaf physiological parameters

[leaf] = LeafPhysiologyParams (params, physcon, leaf);

% --- Plot type: 1 = light. 2 = CO2. 3 = temperature. 4 = vapor pressure deficit

plot_type = 1;

% --- Default conditions

co2conc = 380;                   % Atmospheric CO2 (umol/mol)
relhum = 80;                     % Air relative humidity (%)
air_temp = physcon.tfrz + 25;    % Air temperature (K)
par_sat = 2000;                  % Incident PAR (umol photon/m2/s)
leaf_temp = physcon.tfrz + 25;   % Leaf temperature (K)

% Atmospheric CO2 (umol/mol) and O2 (mmol/mol)

atmos.co2air = co2conc;
atmos.o2air = 0.209 * 1000;

% Leaf absorbed PAR (umol photon/m2 leaf/s)

atmos.par = par_sat;
flux.apar = atmos.par * (1 - leaf.rho(params.vis) - leaf.tau(params.vis));

% Leaf temperature (K) and saturation vapor pressure (Pa)

flux.tleaf = leaf_temp;
[esat_tleaf, desat_tleaf] = satvap ((flux.tleaf-physcon.tfrz));

% Air temperature (K) and vapor pressure (Pa)

atmos.tair = air_temp;
atmos.relhum = relhum;
[esat_tair, desat_tair] = satvap ((atmos.tair-physcon.tfrz));
atmos.eair = esat_tair * (atmos.relhum / 100);
vpd_tleaf = esat_tleaf - atmos.eair;

% Wind (m/s)
% u = 0.01_r8  ! Still air
% u = 0.1_r8   ! Calm - smoke rises vertically
% u = 1.0_r8   ! Light air - smoke drift indicates wind direction
% u = 2.5_r8   ! Light breeze - wind felt on skin, leaves rustle
% u = 5.0_r8   ! Gentle breeze - leaves constantly moving and light flag extended
% u = 10.0_r8  ! Moderate breeze

atmos.wind = 5;

% Atmospheric pressure (Pa) and molar density (mol/m3)

atmos.patm = 101325;
atmos.rhomol = atmos.patm / (physcon.rgas * atmos.tair);

% --- Leaf boundary layer conductances

[flux] = LeafBoundaryLayer (physcon, atmos, leaf, flux);

% ---  Light response at standard CO2, leaf temperature, and VPD

if (plot_type == 1)

   p = 0;
   for i = 0: 10: 2000

      % Set value for PAR

      atmos.par = i;
      flux.apar = atmos.par * (1 - leaf.rho(params.vis) - leaf.tau(params.vis));

      % Calculate photosynthesis and stomatal conductance

      if (leaf.gstyp <= 1)
         [flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);
      elseif (leaf.gstyp == 2)
         [flux] = StomataOptimization (physcon, atmos, leaf, flux);
      end

      % Save data for output

      p = p + 1;
      x1(p) = atmos.par;
      x2(p) = flux.tleaf - physcon.tfrz;
      x3(p) = atmos.co2air;
      x4(p) = (esat_tleaf - atmos.eair) * 0.001;
      x5(p) = flux.hs;
      x6(p) = flux.ci / atmos.co2air;
      x7(p) = flux.ac - flux.rd;
      x8(p) = flux.aj - flux.rd;
      x9(p) = flux.ap - flux.rd;
      x10(p) = flux.an;
      x11(p) = flux.gs;
      x12(p) = flux.gbh;

   end

   % Plot data

   plot(x1,x11)
   xlabel('PAR (\mumol m^{-2} s^{-1})')
   ylabel('Stomatal conductance (mol H_2O m^{-2} s^{-1})')

   % Output data and then clear variables from memory

   A = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12];
   filename = 'light_response.txt';
   fileID = fopen(filename,'w');
   fprintf(fileID,'%6.1f %6.1f %6.1f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n', A);
   fclose(fileID);
   clear x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 A;

   % Reset to default value

   atmos.par = par_sat;
   flux.apar = atmos.par * (1 - leaf.rho(params.vis) - leaf.tau(params.vis));

end

% --- CO2 response at saturated light and standard leaf temperature and VPD

if (plot_type == 2)

   p = 0;
   for i = 100: 10: 1000

      % Set value for CO2

      atmos.co2air = i;

      % Calculate photosynthesis and stomatal conductance

      if (leaf.gstyp <= 1)
         [flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);
      elseif (leaf.gstyp == 2)
         [flux] = StomataOptimization (physcon, atmos, leaf, flux);
      end

      % Save data for output

      p = p + 1;
      x1(p) = atmos.par;
      x2(p) = flux.tleaf - physcon.tfrz;
      x3(p) = atmos.co2air;
      x4(p) = (esat_tleaf - atmos.eair) * 0.001;
      x5(p) = flux.hs;
      x6(p) = flux.ci / atmos.co2air;
      x7(p) = flux.ac - flux.rd;
      x8(p) = flux.aj - flux.rd;
      x9(p) = flux.ap - flux.rd;
      x10(p) = flux.an;
      x11(p) = flux.gs;
      x12(p) = flux.gbh;

   end

   % Plot data

   plot(x3,x11)
   xlabel('c_a (\mumol mol^{-1})')
   ylabel('Stomatal conductance (mol H_2O m^{-2} s^{-1})')

   % Output data and then clear variables from memory

   A = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12];
   filename = 'co2_response.txt';
   fileID = fopen(filename,'w');
   fprintf(fileID,'%6.1f %6.1f %6.1f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n', A);
   fclose(fileID);
   clear x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 A;

   % Reset to default value

   atmos.co2air = co2conc;

end

% --- Temperature response at saturated light and standard CO2 and constant RH or VPD

if (plot_type == 3)

   p = 0;
   for i = 10: 1: 35

      % Set value for leaf temperature and adjust vapor pressure of air so that
      % RH or VPD is constant

      flux.tleaf = physcon.tfrz + i;
      [esat_current, desat_current] = satvap ((flux.tleaf-physcon.tfrz));

      if (leaf.gstyp == 1)
         % Adjust vapor pressure of air so that RH is constant
         atmos.eair = (relhum  / 100) * esat_current;
      else
         % Adjust vapor pressure of air so that VPD is constant
         atmos.eair = esat_current - vpd_tleaf;
      end

      % Change iota based on temperature

      %  if (leaf.gstyp == 2)
      %     ft = @(tl, ha) exp(ha/(physcon.rgas*(physcon.tfrz+25)) * (1-(physcon.tfrz+25)/tl));
      %     cp_tleaf = leaf.cp25 * ft(flux.tleaf, leaf.cpha);
      %     leaf.iota = (leaf.iota25 / 1022) * 3 * cp_tleaf / (1.6 * 2.8 * 2.8) * 100;
      %  end

      % Calculate photosynthesis and stomatal conductance

      if (leaf.gstyp <= 1)
         [flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);
      elseif (leaf.gstyp == 2)
         [flux] = StomataOptimization (physcon, atmos, leaf, flux);
      end

      % Save data for output

      p = p + 1;
      x1(p) = atmos.par;
      x2(p) = flux.tleaf - physcon.tfrz;
      x3(p) = atmos.co2air;
      x4(p) = (esat_current - atmos.eair) * 0.001;
      x5(p) = flux.hs;
      x6(p) = flux.ci / atmos.co2air;
      x7(p) = flux.ac - flux.rd;
      x8(p) = flux.aj - flux.rd;
      x9(p) = flux.ap - flux.rd;
      x10(p) = flux.an;
      x11(p) = flux.gs;
      x12(p) = flux.gbh;

   end

   % Plot data

   plot(x2,x11)
   xlabel('Temperature (^{o}C)')
   ylabel('Stomatal conductance (mol H_2O m^{-2} s^{-1})')

   % Output data and then clear variables from memory

   A = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12];
   filename = 'temperature_response.txt';
   fileID = fopen(filename,'w');
   fprintf(fileID,'%6.1f %6.1f %6.1f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n', A);
   fclose(fileID);
   clear x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 A;

   % Reset to default value

   flux.tleaf = leaf_temp;
   atmos.eair = (relhum  / 100) * esat_tleaf;
   if (leaf.gstyp == 2)
      leaf.iota = leaf.iota25;
   end

end

% --- Vapor pressure response at saturated light and standard CO2 and temperature

if (plot_type == 4)

   p = 0;
   for i = 10: 2: 100

      % Set value for vapor pressure of air

      [esat, desat] = satvap ((atmos.tair-physcon.tfrz));
      atmos.relhum = i;
      atmos.eair = esat * (atmos.relhum / 100);

      % Calculate photosynthesis and stomatal conductance

      if (leaf.gstyp <= 1)
         [flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);
      elseif (leaf.gstyp == 2)
         [flux] = StomataOptimization (physcon, atmos, leaf, flux);
      end

      % Save data for output

      p = p + 1;
      x1(p) = atmos.par;
      x2(p) = flux.tleaf - physcon.tfrz;
      x3(p) = atmos.co2air;
      x4(p) = (esat_tleaf - atmos.eair) * 0.001;
      x5(p) = flux.hs;
      x6(p) = flux.ci / atmos.co2air;
      x7(p) = flux.ac - flux.rd;
      x8(p) = flux.aj - flux.rd;
      x9(p) = flux.ap - flux.rd;
      x10(p) = flux.an;
      x11(p) = flux.gs;
      x12(p) = flux.gbh;

   end

   % Plot data

   plot(x4,x11)
   xlabel('VPD (kPa)')
   ylabel('Stomatal conductance (mol H_2O m^{-2} s^{-1})')

   % Output data and then clear variables from memory

   A = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12];
   filename = 'vpd_response.txt';
   fileID = fopen(filename,'w');
   fprintf(fileID,'%6.1f %6.1f %6.1f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n', A);
   fclose(fileID);
   clear x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 A;

   % Reset to default value

   atmos.relhum = relhum;
   [esat, desat] = satvap ((atmos.tair-physcon.tfrz));
   atmos.eair = esat * atmos.relhum / 100;

end
