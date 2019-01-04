function [flux] = GoudriaanRadiation (omega, Kb, Kbm, Kdm, albcanb, albcand, params, canopy, atmos, flux)

% Compute solar radiation transfer through canopy using Goudriaan parameterization.
% Radiative transfer of Goudriaan (1977), as described by Goudriaan and van Laar
% (1994) and implemented in the plant canopy models of Spitters (1986),
% de Pury and Farquhar (1997), Wang and Leuning (1998), and  Wang (2003).
% This uses the multilayer form of the equations, with constant optical
% properties through the canopy.

% -----------------------------------------------------------------------
% Input
% omega           ! Leaf scattering coefficient for canopy
% Kb              ! Direct beam extinction coefficient for canopy
% Kbm             ! Direct beam extinction coefficient for canopy adjusted for scattering
% Kdm             ! Diffuse extinction coefficient for canopy adjusted for scattering
% albcanb         ! Direct beam albedo above canopy
% albcand         ! Diffuse albedo above canopy
% params.numrad   ! Number of wavebands
% params.npts     ! Number of grid points to process
% params.sun      ! Index for sunlit leaf
% params.sha      ! Index for shaded leaf
% canopy.ntop     ! Index for top leaf layer
% canopy.nbot     ! Index for bottom leaf layer
% canopy.dlai     ! Layer leaf area index (m2/m2)
% canopy.lai      ! Leaf area index of canopy (m2/m2)
% canopy.sumlai   ! Cumulative leaf area index for canopy layer (m2/m2)
% canopy.clumpfac ! Leaf clumping index
% atmos.swskyb    ! Atmospheric direct beam solar radiation (W/m2)
% atmos.swskyd    ! Atmospheric diffuse solar radiation (W/m2)
% flux.fracsun    ! Sunlit fraction of canopy layer
% flux.fracsha    ! Shaded fraction of canopy layer
%
% Output
% flux.swleaf     ! Leaf absorbed solar radiation (W/m2 leaf)
% flux.swveg      ! Absorbed solar radiation, vegetation (W/m2)
% flux.swvegsun   ! Absorbed solar radiation, sunlit canopy (W/m2)
% flux.swvegsha   ! Absorbed solar radiation, shaded canopy (W/m2)
% flux.swsoi      ! Absorbed solar radiation, ground (W/m2)
% flux.albcan     ! Albedo above canopy
% -----------------------------------------------------------------------

% --- Process each waveband (ib) for each grid point (p)

for ib = 1:params.numrad
   for p = 1:params.npts

      % Zero terms that are summed over all layers

      icsun = 0;
      icsha = 0;
      icshad = 0;
      icshabs = 0;
      icsund = 0;
      icsunbs = 0;
      icsunb = 0;

      % Process each canopy layer (iv)

      for iv = canopy.nbot(p):canopy.ntop(p)

         % --- Calculate leaf fluxes. Fluxes are per unit leaf area (W/m2 leaf)
         
         % ild - absorbed diffuse flux per unit leaf area at cumulative LAI, 
         % average for all leaves (J / m2 leaf / s)

         ild = (1 - albcand(p,ib)) * atmos.swskyd(p,ib) * Kdm(p,ib) * canopy.clumpfac(p) ...
             * exp(-Kdm(p,ib) * canopy.sumlai(p,iv) * canopy.clumpfac(p));

         % ilb - absorbed direct beam flux (total with scattering) per unit leaf area 
         % at cumulative LAI, average for all leaves (J / m2 leaf / s)

         ilb = (1 - albcanb(p,ib)) * atmos.swskyb(p,ib) * Kbm(p,ib) * canopy.clumpfac(p) ...
             * exp(-Kbm(p,ib) * canopy.sumlai(p,iv) * canopy.clumpfac(p));

         % ilbb - absorbed direct beam flux (unscattered direct component) per unit leaf area 
         % at cumulative LAI, average for all leaves (J / m2 leaf / s)

         ilbb = (1 - omega(p,ib)) * atmos.swskyb(p,ib) * Kb(p) * canopy.clumpfac(p) ...
              * exp(-Kb(p) * canopy.sumlai(p,iv) * canopy.clumpfac(p));

         % ilbs - absorbed direct beam flux (scattered direct component) per unit leaf area 
         % at cumulative LAI, average for all leaves (J / m2 leaf / s)

         ilbs = ilb - ilbb;

         % ilsha - total absorbed flux (shaded leaves) per unit shaded leaf area 
         % at cumulative LAI (J / m2 leaf / s)

         ilsha = ild + ilbs;

         % ilsun - total absorbed flux (sunlit leaves) per unit sunlit leaf area 
         % at cumulative LAI (J / m2 leaf / s)

         ilsun = ilsha + Kb(p) * (1 - omega(p,ib)) * atmos.swskyb(p,ib);

         % Save solar radiation absorbed by sunlit and shaded leaves

         flux.swleaf(p,iv,params.sun,ib) = ilsun;
         flux.swleaf(p,iv,params.sha,ib) = ilsha;

         % --- Canopy summation and soil absorption. Fluxes are per unit ground area (W/m2 ground area)

         % icsun - absorbed solar radiation, sunlit canopy (W/m2)
         % icsha - absorbed solar radiation, shaded canopy (W/m2)

         icsun = icsun + ilsun * flux.fracsun(p,iv) * canopy.dlai(p,iv);
         icsha = icsha + ilsha * flux.fracsha(p,iv) * canopy.dlai(p,iv);

         % icshad  - diffuse radiation absorbed by shaded leaves (W/m2)
         % icshabs - scattered direct beam radiation absorbed by shaded leaves (W/m2)

         icshad = icshad + ild * flux.fracsha(p,iv) * canopy.dlai(p,iv);
         icshabs = icshabs + ilbs * flux.fracsha(p,iv) * canopy.dlai(p,iv);

         % icsund  - diffuse radiation absorbed by sunlit leaves (W/m2)
         % icsunbs - scattered direct beam radiation absorbed by sunlit leaves (W/m2)
         % icsunb  - direct beam radiation absorbed by sunlit leaves (W/m2)

         icsund = icsund + ild * flux.fracsun(p,iv) * canopy.dlai(p,iv);
         icsunbs = icsunbs + ilbs * flux.fracsun(p,iv) * canopy.dlai(p,iv);
         icsunb = icsunb + Kb(p) * (1 - omega(p,ib)) * atmos.swskyb(p,ib) * flux.fracsun(p,iv) * canopy.dlai(p,iv);

      end

      % Solar radiation absorbed by vegetation (W/m2)

      sabv = icsun + icsha;

      % Solar radiation absorbed by ground (W/m2)

      sabg = atmos.swskyb(p,ib) * (1 - albcanb(p,ib)) * exp(-Kbm(p,ib)*canopy.lai(p)*canopy.clumpfac(p)) + ...
             atmos.swskyd(p,ib) * (1 - albcand(p,ib)) * exp(-Kdm(p,ib)*canopy.lai(p)*canopy.clumpfac(p));

      % Conservation check: absorbed = incoming - outgoing
      % This is not valid, because the numerical integration of leaf fluxes does not equal the analytical
      % solution (unless dlai is very small)

      suminc = atmos.swskyb(p,ib) + atmos.swskyd(p,ib);
      sumabs = sabv + sabg;
      sumref = atmos.swskyb(p,ib) * albcanb(p,ib) + atmos.swskyd(p,ib) * albcand(p,ib);

      err = suminc - (sumref + sumabs);
      err = 0;
      if (abs(err) > 1e-03)
         fprintf('err = %15.5f\n',err)
         fprintf('suminc = %15.5f\n',suminc)
         fprintf('sumref = %15.5f\n',sumref)
         fprintf('sumabs = %15.5f\n',sumabs)
         error ('GoudriaanRadiation: Solar radiation conservation error')
      end

      % --- Save necessary radiative fluxes

      % Albedo

      suminc = atmos.swskyb(p,ib) + atmos.swskyd(p,ib);
      sumref = atmos.swskyb(p,ib) * albcanb(p,ib) + atmos.swskyd(p,ib) * albcand(p,ib);
      if (suminc > 0)
         flux.albcan(p,ib) = sumref / suminc;
      else
         flux.albcan(p,ib) = 0;
      end

      % Solar radiation absorbed by canopy

      flux.swveg(p,ib) = sabv;
      flux.swvegsun(p,ib) = icsun;
      flux.swvegsha(p,ib) = icsha;

      % Solar radiation absorbed by ground (soil)

      flux.swsoi(p,ib) = sabg;

   end
end
