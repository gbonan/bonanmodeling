% Supplemental program 11.1

% --------------------------------------------------------------------
% Calculate C3 photosynthesis in relation to PAR, CO2, and temperature
% --------------------------------------------------------------------

% --- Waveband indices for visible and near-infrared

params.vis = 1; params.nir = 2;

% --- Physical constants

physcon.tfrz = 273.15;                % Freezing point of water (K)
physcon.rgas = 8.31446;               % Universal gas constant (J/K/mol)

% --- Set leaf physiology variables

% Photosynthesis co-limitation: 0 = no. 1 = yes

  leaf.colim = 1;

% Leaf physiological parameters

[leaf] = LeafPhysiologyParams (params, physcon, leaf);

% --- Plot type: 1 = light. 2 = CO2. 3 = temperature

plot_type = 1;

% --- Default conditions

co2conc = 380;                   % Atmospheric CO2 (umol/mol)
par_sat = 2000;                  % Incident PAR (umol photon/m2/s)
leaf_temp = physcon.tfrz + 25;   % Leaf temperature (K)

% Atmospheric CO2 (umol/mol) and O2 (mmol/mol)

atmos.co2air = co2conc;
atmos.o2air = 0.209 * 1000;

% Leaf absorbed PAR (umol photon/m2 leaf/s)

atmos.par = par_sat;
flux.apar = atmos.par * (1 - leaf.rho(params.vis) - leaf.tau(params.vis));

% Leaf temperature (K)

flux.tleaf = leaf_temp;

% ---  Light response at standard CO2 and leaf temperature

if (plot_type == 1)

   p = 0;
   for i = 0: 10: 2000

      % Set value for PAR

      atmos.par = i;
      flux.apar = atmos.par * (1 - leaf.rho(params.vis) - leaf.tau(params.vis));

      % Calculate photosynthesis

      [flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);

      % Save data for output

      p = p + 1;
      x1(p) = atmos.par;
      x2(p) = flux.tleaf - physcon.tfrz;
      x3(p) = atmos.co2air;
      x4(p) = flux.ci / atmos.co2air;
      x5(p) = flux.ac - flux.rd;
      x6(p) = flux.aj - flux.rd;
      x7(p) = flux.an;

   end

   % Plot data

   plot(x1,x5,'b',x1,x6,'r',x1,x7,'g')
   title('C_3 photosynthesis')
   xlabel('PAR (\mumol m^{-2} s^{-1})')
   ylabel('Photosynthesis (\mumol CO_2 m^{-2} s^{-1})')
   legend('A_c','A_j','A_n','Location','best')

   % Output data and then clear variables from memory

   A = [x1; x2; x3; x4; x5; x6; x7];
   filename = 'light_response.txt';
   fileID = fopen(filename,'w');
   fprintf(fileID,'%9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n', A);
   fclose(fileID);
   clear x1 x2 x3 x4 x5 x6 x7 A;

   % Reset to default value

   atmos.par = par_sat;
   flux.apar = atmos.par * (1 - leaf.rho(params.vis) - leaf.tau(params.vis));

end

% --- CO2 response at saturated light and standard leaf temperature

if (plot_type == 2)

   p = 0;
   for i = 40: 10: 1000

      % Set value for CO2

      atmos.co2air = i / 0.7;

      % Calculate photosynthesis

      [flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);

      % Save data for output

      p = p + 1;
      x1(p) = atmos.par;
      x2(p) = flux.tleaf - physcon.tfrz;
      x3(p) = atmos.co2air;
      x4(p) = flux.ci / atmos.co2air;
      x5(p) = flux.ac - flux.rd;
      x6(p) = flux.aj - flux.rd;
      x7(p) = flux.an;

   end

   % Plot data

   plot(x3,x5,'b',x3,x6,'r',x3,x7,'g')
   title('C_3 photosynthesis')
   xlabel('c_i (\mumol mol^{-1})')
   ylabel('Photosynthesis (\mumol CO_2 m^{-2} s^{-1})')
   legend('A_c','A_j','A_n','Location','best')

   % Output data and then clear variables from memory

   A = [x1; x2; x3; x4; x5; x6; x7];
   filename = 'co2_response.txt';
   fileID = fopen(filename,'w');
   fprintf(fileID,'%9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n', A);
   fclose(fileID);
   clear x1 x2 x3 x4 x5 x6 x7 A;

   % Reset to default value

   atmos.co2air = co2conc;

end

% --- Temperature response at saturated light and standard CO2

if (plot_type == 3)

   p = 0;
   for i = 10: 1: 35

      % Set value for leaf temperature

      flux.tleaf = physcon.tfrz + i;

      % Calculate photosynthesis

      [flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux);

      % Save data for output

      p = p + 1;
      x1(p) = atmos.par;
      x2(p) = flux.tleaf - physcon.tfrz;
      x3(p) = atmos.co2air;
      x4(p) = flux.ci / atmos.co2air;
      x5(p) = flux.ac - flux.rd;
      x6(p) = flux.aj - flux.rd;
      x7(p) = flux.an;

   end

   % Plot data

   plot(x2,x5,'b',x2,x6,'r',x2,x7,'g')
   title('C_3 photosynthesis')
   xlabel('Temperature (^{o}C)')
   ylabel('Photosynthesis (\mumol CO_2 m^{-2} s^{-1})')
   legend('A_c','A_j','A_n','Location','best')

   % Output data and then clear variables from memory

   A = [x1; x2; x3; x4; x5; x6; x7];
   filename = 'temperature_response.txt';
   fileID = fopen(filename,'w');
   fprintf(fileID,'%9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n', A);
   fclose(fileID);
   clear x1 x2 x3 x4 x5 x6 x7 A;

   % Reset to default value

   flux.tleaf = leaf_temp;

end
