function [flux] = TwoStreamRadiation (omega, avmu, betad, betab, Kb, params, canopy, atmos, flux)

% Compute solar radiation transfer through canopy using the two-stream approximation

% -----------------------------------------------------------------------
% Input
% omega           ! Leaf scattering coefficient for canopy
% avmu            ! Average inverse diffuse optical depth per unit leaf area
% betad           ! Upscatter parameter for diffuse radiation
% betab           ! Upscatter parameter for direct beam radiation
% Kb              ! Optical depth of direct beam per unit leaf area (direct beam extinction coefficient for canopy)
% params.numrad   ! Number of wavebands
% params.npts     ! Number of grid points to process
% params.sun      ! Index for sunlit leaf
% params.sha      ! Index for shaded leaf
% canopy.ntop     ! Index for top leaf layer
% canopy.nbot     ! Index for bottom leaf layer
% canopy.lai      ! Leaf area index of canopy (m2/m2)
% canopy.sumlai   ! Cumulative leaf area index for canopy layer (m2/m2)
% canopy.dlai     ! Layer leaf area index (m2/m2)
% canopy.clumpfac ! Leaf clumping index
% atmos.swskyb    ! Atmospheric direct beam solar radiation (W/m2)
% atmos.swskyd    ! Atmospheric diffuse solar radiation (W/m2)
% flux.fracsun    ! Sunlit fraction of canopy layer
% flux.fracsha    ! Shaded fraction of canopy layer
% flux.albsoib    ! Direct beam albedo of ground (soil)
% flux.albsoid    ! Diffuse albedo of ground (soil)
%
% Output
% flux.swleaf     ! Leaf absorbed solar radiation (W/m2 leaf)
% flux.swveg      ! Absorbed solar radiation, vegetation (W/m2)
% flux.swvegsun   ! Absorbed solar radiation, sunlit canopy (W/m2)
% flux.swvegsha   ! Absorbed solar radiation, shaded canopy (W/m2)
% flux.swsoi      ! Absorbed solar radiation, ground (W/m2)
% flux.albcan     ! Albedo above canopy
% -----------------------------------------------------------------------

% --- Process each waveband for each grid point

for ib = 1:params.numrad
   for p = 1:params.npts

      % --- Canopy fluxes using total canopy lai

      % Common terms

      b = (1 - (1 - betad(p,ib)) * omega(p,ib)) / avmu(p);
      c = betad(p,ib) * omega(p,ib) / avmu(p);
      h = sqrt(b*b - c*c);
      u = (h - b - c) / (2 * h);
      v = (h + b + c) / (2 * h);
      d = omega(p,ib) * Kb(p) * atmos.swskyb(p,ib) / (h*h - Kb(p)*Kb(p));
      g1 = (betab(p,ib) * Kb(p) - b * betab(p,ib) - c * (1 - betab(p,ib))) * d;
      g2 = ((1 - betab(p,ib)) * Kb(p) + c * betab(p,ib) + b * (1 - betab(p,ib))) * d;
      s1 = exp(-h * canopy.lai(p) * canopy.clumpfac(p));
      s2 = exp(-Kb(p) * canopy.lai(p) * canopy.clumpfac(p));

      % Direct beam radiation

      num1 = v * (g1 + g2 * flux.albsoid(p,ib) + flux.albsoib(p,ib) * atmos.swskyb(p,ib)) * s2;
      num2 = g2 * (u + v * flux.albsoid(p,ib)) * s1;
      den1 = v * (v + u * flux.albsoid(p,ib)) / s1;
      den2 = u * (u + v * flux.albsoid(p,ib)) * s1;
      n2b = (num1 - num2) / (den1 - den2);
      n1b = (g2 - n2b * u) / v;

      a1b = -g1 *      (1 - s2*s2) / (2 * Kb(p)) + ...
             n1b * u * (1 - s2*s1) / (Kb(p) + h) + n2b * v * (1 - s2/s1) / (Kb(p) - h);
      a2b =  g2 *      (1 - s2*s2) / (2 * Kb(p)) - ...
             n1b * v * (1 - s2*s1) / (Kb(p) + h) - n2b * u * (1 - s2/s1) / (Kb(p) - h);

      % iupwb0    - Direct beam flux scattered upward (reflected) above canopy (W/m2)
      % iupwb     - Direct beam flux scattered upward at the canopy depth (W/m2)
      % idwnb     - Direct beam flux scattered downward below canopy (W/m2)
      % iabsb     - Direct beam flux absorbed by canopy (W/m2)
      % iabsb_sun - Direct beam flux absorbed by sunlit canopy (W/m2)
      % iabsb_sha - Direct beam flux absorbed by shaded canopy (W/m2)

      iupwb0 = -g1 + n1b * u + n2b * v;
      iupwb = -g1 * s2 + n1b * u * s1 + n2b * v / s1;
      idwnb = g2 * s2 - n1b * v * s1 - n2b * u / s1;
      iabsb = atmos.swskyb(p,ib) * (1 - s2) - iupwb0 + iupwb - idwnb;
      iabsb_sun = (1 - omega(p,ib)) ...
         * ((1 - s2) * atmos.swskyb(p,ib) + 1 / avmu(p) * (a1b + a2b) * canopy.clumpfac(p));
      iabsb_sha = iabsb - iabsb_sun;

      % Diffuse radiation
 
      num = atmos.swskyd(p,ib) * (u + v * flux.albsoid(p,ib)) * s1;
      den1 = v * (v + u * flux.albsoid(p,ib)) / s1;
      den2 = u * (u + v * flux.albsoid(p,ib)) * s1;
      n2d = num / (den1 - den2);
      n1d = -(atmos.swskyd(p,ib) + n2d * u) / v;

      a1d =  n1d * u * (1 - s2*s1) / (Kb(p) + h) + n2d * v * (1 - s2/s1) / (Kb(p) - h);
      a2d = -n1d * v * (1 - s2*s1) / (Kb(p) + h) - n2d * u * (1 - s2/s1) / (Kb(p) - h);

      % iupwd0    - Diffuse flux scattered upward (reflected) above canopy (W/m2)
      % iupwd     - Diffuse flux scattered upward at the canopy depth (W/m2)
      % idwnd     - Diffuse flux scattered downward below canopy (W/m2)
      % iabsd     - Diffuse flux absorbed by canopy (W/m2)
      % iabsd_sun - Diffuse flux absorbed by sunlit canopy (W/m2)
      % iabsd_sha - Diffuse flux absorbed by shaded canopy (W/m2)

      iupwd0 = n1d * u + n2d * v;
      iupwd = n1d * u * s1 + n2d * v / s1;
      idwnd = -n1d * v * s1 - n2d * u / s1;
      iabsd = atmos.swskyd(p,ib) - iupwd0 + iupwd - idwnd;
      iabsd_sun = (1 - omega(p,ib)) / avmu(p) * (a1d + a2d) * canopy.clumpfac(p);
      iabsd_sha = iabsd - iabsd_sun;

      % --- Save necessary radiative fluxes

      % Albedo

      suminc = atmos.swskyb(p,ib) + atmos.swskyd(p,ib);
      sumref = iupwb0 + iupwd0;
      if (suminc > 0)
         flux.albcan(p,ib) = sumref / suminc;
      else
         flux.albcan(p,ib) = 0;
      end

      % Solar radiation absorbed by canopy

      flux.swveg(p,ib) = iabsb +  iabsd;
      flux.swvegsun(p,ib) = iabsb_sun + iabsd_sun;
      flux.swvegsha(p,ib) = iabsb_sha + iabsd_sha;

      % Solar radiation absorbed by ground (soil)

      dir = atmos.swskyb(p,ib) * s2 * (1 - flux.albsoib(p,ib));
      dif = (idwnb + idwnd) * (1 - flux.albsoid(p,ib));
      flux.swsoi(p,ib) = dir + dif;

      % --- Conservation check: total incident = total reflected + total absorbed

      suminc = atmos.swskyb(p,ib) + atmos.swskyd(p,ib);
      sumref = iupwb0 + iupwd0;
      sumabs = flux.swveg(p,ib) + flux.swsoi(p,ib);

      err = suminc - (sumabs + sumref);
      if (abs(err) > 1e-06)
         fprintf('suminc = %15.5f\n',suminc)
         fprintf('sumref = %15.5f\n',sumref)
         fprintf('sumabs = %15.5f\n',sumabs)
         error ('TwoStreamRadiation: Total solar radiation conservation error')
      end

      % --- Repeat two-stream calculations for each leaf layer to calculate leaf fluxes

      icsun(ib) = 0;
      icsha(ib) = 0;

      for iv = canopy.nbot(p):canopy.ntop(p)

         % s1 and s2 depend on cumulative lai

         s1 = exp(-h * canopy.sumlai(p,iv) * canopy.clumpfac(p));
         s2 = exp(-Kb(p) * canopy.sumlai(p,iv) * canopy.clumpfac(p));

         % ilbb - absorbed direct beam flux (unscattered direct component) per unit leaf area
         % at cumulative LAI, average for all leaves (J / m2 leaf / s)

         ilbb = (1 - omega(p,ib)) * Kb(p) * atmos.swskyb(p,ib) * s2;

         % ilbs - absorbed direct beam flux (scattered direct component) per unit leaf area
         % at cumulative LAI, average for all leaves (J / m2 leaf / s)

         diupwb = Kb(p) * g1 * s2 - h * n1b * u * s1 + h * n2b * v / s1;
         didwnb = -Kb(p) * g2 * s2 + h * n1b * v * s1 - h * n2b * u / s1;
         ilbs = (omega(p,ib) * Kb(p) * atmos.swskyb(p,ib) * s2 + (diupwb - didwnb)) * canopy.clumpfac(p);

         % ild - absorbed diffuse flux per unit leaf area at cumulative LAI,
         % average for all leaves (J / m2 leaf / s)

         diupwd = -h * n1d * u * s1 + h * n2d * v / s1;
         didwnd = h * n1d * v * s1 - h * n2d * u / s1;
         ild = (diupwd - didwnd) * canopy.clumpfac(p);

         % Save leaf fluxes per unit sunlit and shaded leaf area (W/m2 leaf)

         flux.swleaf(p,iv,params.sun,ib) = ilbb / flux.fracsun(p,iv) + (ilbs + ild);
         flux.swleaf(p,iv,params.sha,ib) = ilbs + ild;

         icsun(ib) = icsun(ib) + (ilbb + (ilbs + ild)*flux.fracsun(p,iv)) * canopy.dlai(p,iv);
         icsha(ib) = icsha(ib) + (ilbs + ild)*flux.fracsha(p,iv) * canopy.dlai(p,iv);

      end   % end canopy loop
   end      % end grid point loop
end         % end waveband loop

% --- Adjust leaf fluxes as needed. The sum of the fluxes for sunlit and shaded
% leaves should equal the total absorbed by the canopy, but may not because of
% inaccuracies in the flux derivatives (this is a small error if the dlai increment
% is small). Normalize these fluxes to sum to the canopy absorption.

for ib = 1:params.numrad
   for p = 1:params.npts

      % Sum canopy absorption (W/m2 ground) using leaf fluxes per unit sunlit
      % and shaded leaf area (W/m2 leaf)

      sumabs = 0;
      sumabs_sun = 0;
      sumabs_sha = 0;
      for iv = canopy.nbot(p):canopy.ntop(p)
         sun = flux.swleaf(p,iv,params.sun,ib) * flux.fracsun(p,iv) * canopy.dlai(p,iv);
         sha = flux.swleaf(p,iv,params.sha,ib) * flux.fracsha(p,iv) * canopy.dlai(p,iv);
         sumabs = sumabs + sun + sha;
         sumabs_sun = sumabs_sun + sun;
         sumabs_sha = sumabs_sha + sha;
      end

      % Normalize profile

      if (sumabs > 0)
         for iv = canopy.nbot(p):canopy.ntop(p)
            flux.swleaf(p,iv,params.sun,ib) = flux.swleaf(p,iv,params.sun,ib) * flux.swveg(p,ib) / sumabs;
            flux.swleaf(p,iv,params.sha,ib) = flux.swleaf(p,iv,params.sha,ib) * flux.swveg(p,ib) / sumabs;
         end
      end

   end

end
