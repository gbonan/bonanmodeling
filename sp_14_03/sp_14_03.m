% Supplemental program 14.3

% --------------------------------------------
% Calculate and graph light profiles in canopy
% --------------------------------------------

% --- Model parameters

params.numrad = 2;         % Number of wavebands (visible, near-infrared)
params.vis = 1;            % Array index for visible waveband
params.nir = 2;            % Array index for near-infrared waveband
params.sun = 1;            % Array index for sunlit leaf
params.sha = 2;            % Array index for shaded leaf
params.npts = 1;           % Number of grid points to process

% --- Model options

  light = 'Norman';        % Use Norman radiative transfer
% light = 'Goudriaan';     % Use Goudriaan radiative transfer
% light = 'TwoStream';     % Use two-stream approximation radiative transfer

  canopy_type = 'dense';   % High leaf area index
% canopy_type = 'sparse';  % Low leaf area index

% --- Define plant canopy

for p = 1:params.npts

   % Set canopy LAI, layer LAI increment, and number of layers

   switch canopy_type
   case 'dense'
      lai_inc = 0.1;                                % Leaf area index for each layer
      canopy.lai(p) = 6;                            % Leaf area index of canopy (m2/m2)
   case 'sparse'
      lai_inc = 0.05;                               % Leaf area index for each layer
      canopy.lai(p) = 1;                            % Leaf area index of canopy (m2/m2)
   end
   canopy.nveg(p) = round(canopy.lai(p) / lai_inc); % Number of leaf layers in canopy

   % Minimum number of layers for Norman radiation

   switch light
   case 'Norman'
      if (canopy.nveg(p) < 10)
         canopy.nveg(p) = 10;
         lai_inc = canopy.lai(p) / canopy.nveg(p);
      end
   end

   % Set array indices for canopy layers

   canopy.nsoi(p) = 1;                                     % First layer is soil
   canopy.nbot(p) = canopy.nsoi(p) + 1;                    % Bottom leaf layer
   canopy.ntop(p) = canopy.nbot(p) + canopy.nveg(p) - 1;   % Top leaf layer

   % Set LAI of each layer

   for iv = canopy.nbot(p):canopy.ntop(p)
      canopy.dlai(p,iv) = lai_inc;
   end

   % Cumulative leaf area index (from canopy top) at mid-layer

   for iv = canopy.ntop(p): -1: canopy.nbot(p)
      if (iv == canopy.ntop(p))
         canopy.sumlai(p,iv) = 0.5 * canopy.dlai(p,iv);
      else
         canopy.sumlai(p,iv) = canopy.sumlai(p,iv+1) + canopy.dlai(p,iv);
      end
   end

   % Clumping index

   canopy.clumpfac(p) = 1;

end

% --- Atmospheric solar radiation. Solar radiation is given as a unit of visible radiation
% and a unit of near-infrared radiation, both split into direct and diffuse components.

for p = 1:params.npts
   atmos.solar_zenith(p) = 30 * (pi / 180);     % Solar zentih angle (radians)
   atmos.swskyb(p,params.vis) = 0.8;            % Direct beam solar radiation for visible waveband (W/m2)
   atmos.swskyd(p,params.vis) = 0.2;            % Diffuse solar radiation for visible waveband (W/m2)
   atmos.swskyb(p,params.nir) = 0.8;            % Direct beam solar radiation for near-infrared waveband (W/m2)
   atmos.swskyd(p,params.nir) = 0.2;            % Diffuse solar radiation for near-infrared waveband (W/m2)
end

% --- Leaf optical properties

for p = 1:params.npts
   rho(p,params.vis) = 0.10;                    % Leaf reflectance (visible)
   rho(p,params.nir) = 0.45;                    % Leaf reflectance (near-infrared)
   tau(p,params.vis) = 0.05;                    % Leaf transmittance (visible)
   tau(p,params.nir) = 0.25;                    % Leaf transmittance (near-infrared)
   for ib = 1:params.numrad
      omega(p,ib) = rho(p,ib) + tau(p,ib);      % Leaf scattering coefficient for canopy
   end
end

% --- Soil albedo

for p = 1:params.npts
   flux.albsoib(p,params.vis) = 0.1;                         % Direct beam albedo of ground (visible)
   flux.albsoid(p,params.vis) = flux.albsoib(p,params.vis);  % Diffuse albedo of ground (visible)
   flux.albsoib(p,params.nir) = 0.2;                         % Direct beam albedo of ground (near-infrared)
   flux.albsoid(p,params.nir) = flux.albsoib(p,params.nir);  % Diffuse albedo of ground near-infrared)
end

% --- Direct beam extinction coefficient

% xl - departure of leaf angle from spherical orientation

xl = 0.00;

% Kb - direct beam extinction coefficient for canopy

for p = 1:params.npts

   % -0.4 <= xl <= 0.6

   chil(p) = min(max(xl, -0.4), 0.6);

  % Prevent near-zero xl for two-stream radiation

   switch light
   case 'TwoStream'
      if (abs(chil(p)) <= 0.01)
         chil(p) = 0.01;
      end
   end

   % Terms in Ross-Goudriaan function for gdir

   phi1(p) = 0.5 - 0.633 * chil(p) - 0.330 * chil(p)*chil(p);
   phi2(p) = 0.877 * (1 - 2 * phi1(p));

   % Relative projected area of leaf in the direction of solar beam

   cosz = cos(atmos.solar_zenith(p));
   gdir(p) = phi1(p) + phi2(p) * cosz;

   % Direct beam extinction coefficient

   Kb(p) = gdir(p) / cosz;

   % Prevent large Kb at low sun angle

   Kb(p) = min(Kb(p), 20);

end

% --- Sunlit and shaded portions of canopy

for p = 1:params.npts

   % Sunlit and shaded fraction of leaf layer

   for iv = canopy.nbot(p):canopy.ntop(p)
      flux.fracsun(p,iv) = canopy.clumpfac(p) * exp(-Kb(p) * canopy.sumlai(p,iv) * canopy.clumpfac(p));
      flux.fracsha(p,iv) = 1 - flux.fracsun(p,iv);
   end

   % Sunlit and shaded leaf area index for canopy

   laisun(p) = (1 - exp(-Kb(p) * canopy.lai(p) * canopy.clumpfac(p))) / Kb(p);
   laisha(p) = canopy.lai(p) - laisun(p);

end

% --- Unique parameters for Norman radiative transfer

% tb - exponential transmittance of direct beam radiation through a single leaf layer

for p = 1:params.npts
   for iv = canopy.nbot(p):canopy.ntop(p)
      tb(p,iv) = exp(-Kb(p) * canopy.dlai(p,iv) * canopy.clumpfac(p));
   end
end

% td - exponential transmittance of diffuse radiation through a single leaf layer
% with thickness dlai, estimated for nine sky angles in increments of 10 degrees

for p = 1:params.npts
   for iv = canopy.nbot(p):canopy.ntop(p)
      td(p,iv) = 0;

      for j = 1:9

         % Sky angles (5, 15, 25, 35, 45, 55, 65, 75, 85)

         angle = (5 + (j - 1) * 10) * pi / 180;

         % Relative projected area of leaf in the direction of sky angle

         gdirj = phi1(p) + phi2(p) * cos(angle);

         % Sum transmittance

         td(p,iv) = td(p,iv) ...
         + exp(-gdirj / cos(angle) * canopy.dlai(p,iv) * canopy.clumpfac(p)) * sin(angle) * cos(angle);

      end
      td(p,iv) = td(p,iv) * 2 * (10 * pi / 180);

   end
end

% tbcum - direct beam transmittance uses cumulative lai above layer i to
% give unscattered direct beam onto layer i

for p = 1:params.npts
   cumlai = 0;
   iv = canopy.ntop(p);
   tbcum(p,iv) = 1;
   for iv = canopy.ntop(p): -1: canopy.nbot(p)
      cumlai = cumlai + canopy.dlai(p,iv);
      tbcum(p,iv-1) = exp(-Kb(p) * cumlai * canopy.clumpfac(p));
   end
end

% --- Unique parameters for Goudriaan radiative transfer

% Kd - diffuse extinction coefficient for canopy, estimated for nine sky angles
% in increments of 10 degrees

for p = 1:params.npts
   Kd(p) = 0;

   for j = 1:9

      % Sky angles (5, 15, 25, 35, 45, 55, 65, 75, 85)

      angle = (5 + (j - 1) * 10) * pi / 180;

      % Relative projected area of leaf in the direction of sky angle

      gdirj = phi1(p) + phi2(p) * cos(angle);

      % Sum transmittance

      Kd(p) = Kd(p) + exp(-gdirj / cos(angle) * canopy.lai(p) * canopy.clumpfac(p)) * sin(angle) * cos(angle);

   end

   Kd(p) = Kd(p) * 2 * (10 * pi / 180);

   % Convert transmittance to extinction coefficient

   if (canopy.lai(p) > 0)
      Kd(p) = -log(Kd(p)) / (canopy.lai(p) * canopy.clumpfac(p));
   else
      Kd(p) = 0;
   end

end

for ib = 1:params.numrad
   for p = 1:params.npts

      % Adjust Kb and Kd for scattering

      Kbm(p,ib) = Kb(p) * sqrt(1 - omega(p,ib));
      Kdm(p,ib) = Kd(p) * sqrt(1 - omega(p,ib));

      % albvegh - vegetation albedo for horizontal leaves

      albvegh = (1 - sqrt(1 - omega(p,ib))) / (1 + sqrt(1 - omega(p,ib)));

      % albvegb - direct beam vegetation albedo for non-horizontal leaves
      
      albvegb = 2 * Kb(p) / (Kb(p) + Kd(p)) * albvegh;

      % albvegd - diffuse vegetation albedo for non-horizontal leaves, calculated by summing albedo over 9 sky angles

      albvegd = 0;
      for j = 1:9

         % Sky angles (5, 15, 25, 35, 45, 55, 65, 75, 85)

         angle = (5 + (j - 1) * 10) * pi / 180;

         % Relative projected area of leaf in the direction of sky angle

         gdirj = phi1(p) + phi2(p) * cos(angle);

         % Kb for sky angle j

         Kbj = gdirj / cos(angle);

         % Direct beam albedo for sky angle j

         albvegbj = 2 * Kbj / (Kbj + Kd(p)) * albvegh;

         % Sum albedo

         albvegd = albvegd + albvegbj * sin(angle) * cos(angle);

      end
      albvegd = albvegd * 2 * (10 * pi / 180);
   
      % Effective canopy albedo, including soil
      % albcanb - direct beam albedo above canopy
      % albcand - diffuse albedo above canopy
      
      albcanb(p,ib) = albvegb ...
      + (flux.albsoib(p,ib) - albvegb) * exp(-2 * Kbm(p,ib) * canopy.lai(p) * canopy.clumpfac(p));
      albcand(p,ib) = albvegd ...
      + (flux.albsoid(p,ib) - albvegd) * exp(-2 * Kdm(p,ib) * canopy.lai(p) * canopy.clumpfac(p));

   end
end

% --- Two-stream parameters

for p = 1:params.npts

   % avmu - average inverse diffuse optical depth per unit leaf area

   avmu(p) = ( 1 - phi1(p)/phi2(p) * log((phi1(p)+phi2(p))/phi1(p)) ) / phi2(p);

   % Upscatter parameters

   for ib = 1:params.numrad

      % betad - upscatter parameter for diffuse radiation

      betad(p,ib) = 0.5 / omega(p,ib) * ( rho(p,ib) + tau(p,ib) + (rho(p,ib)-tau(p,ib)) * ((1+chil(p))/2)^2 );

      % betab - upscatter parameter for direct beam radiation

      cosz = cos(atmos.solar_zenith(p));
      tmp0 = gdir(p) + phi2(p) * cosz;
      tmp1 = phi1(p) * cosz;
      tmp2 = 1 - tmp1/tmp0 * log((tmp1+tmp0)/tmp1);
      asu = 0.5 * omega(p,ib) * gdir(p) / tmp0 * tmp2;
      betab(p,ib) = (1 + avmu(p)*Kb(p)) / (omega(p,ib)*avmu(p)*Kb(p)) * asu;

   end
end

% --- Light profile through canopy

switch light
   case 'Norman'
   [flux] = NormanRadiation (rho, tau, omega, td, tb, tbcum, params, canopy, atmos, flux);
   case 'Goudriaan'
   [flux] = GoudriaanRadiation (omega, Kb, Kbm, Kdm, albcanb, albcand, params, canopy, atmos, flux);
   case 'TwoStream'
   [flux] = TwoStreamRadiation (omega, avmu, betad, betab, Kb, params, canopy, atmos, flux);
end

% Absorbed PAR per unit sunlit and shaded leaf area (umol photon/m2 leaf/s)

for p = 1:params.npts
   for iv = canopy.nbot(p):canopy.ntop(p)
      flux.apar(p,iv,params.sun) = flux.swleaf(p,iv,params.sun,params.vis) * 4.6;
      flux.apar(p,iv,params.sha) = flux.swleaf(p,iv,params.sha,params.vis) * 4.6;
   end
end

% --- Write formatted output to file, from top layer to bottom layer

p = 1;
for iv = canopy.ntop(p): -1: canopy.nbot(p)
   j = iv - 1;
   k = canopy.nveg(p) - j + 2;
   a1(j) = j;
   a2(j) = canopy.sumlai(p,k);
   a3(j) = flux.swleaf(p,k,params.sun,params.vis);
   a4(j) = flux.swleaf(p,k,params.sha,params.vis);
   a5(j) = flux.swleaf(p,k,params.sun,params.nir);
   a6(j) = flux.swleaf(p,k,params.sha,params.nir);
end

A = [a1; a2; a3; a4; a5; a6];

fileID = fopen('data.txt','w');
fprintf(fileID,'%12s %12s %12s %12s %12s %12s\n','layer','lai','sun','sha','sun','sha');
fprintf(fileID,'%12.0f %12.4f %12.6f %12.6f %12.6f %12.6f\n', A);
fclose(fileID);

% --- Make graph

plot(a3,a2,'b-',a4,a2,'r-',a5,a2,'b--',a6,a2,'r--')
set(gca,'YDir','reverse')
title(light)
xlabel('Fraction absorbed')
ylabel('Cumulative leaf area index')
legend('sun, vis','sha, vis','sun, nir','sha, nir','Location','best')
