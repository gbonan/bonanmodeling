function [flux] = NormanRadiation (rho, tau, omega, td, tb, tbcum, params, canopy, atmos, flux)

% Compute solar radiation transfer through canopy using Norman (1979)

% -----------------------------------------------------------------------
% Input
% rho            ! Leaf reflectance
% tau            ! Leaf transmittance
% omega          ! Leaf scattering coefficient
% td             ! Exponential transmittance of diffuse radiation through a single leaf layer
% tb             ! Exponential transmittance of direct beam radiation through a single leaf layer
% tbcum          ! Cumulative exponential transmittance of direct beam onto a canopy layer
% params.numrad  ! Number of wavebands
% params.npts    ! Number of grid points to process
% params.sun     ! Index for sunlit leaf
% params.sha     ! Index for shaded leaf
% canopy.ntop    ! Index for top leaf layer
% canopy.nbot    ! Index for bottom leaf layer
% canopy.nsoi    ! First canopy layer is soil
% canopy.dlai    ! Layer leaf area index (m2/m2)
% atmos.swskyb   ! Atmospheric direct beam solar radiation (W/m2)
% atmos.swskyd   ! Atmospheric diffuse solar radiation (W/m2)
% flux.fracsun   ! Sunlit fraction of canopy layer
% flux.fracsha   ! Shaded fraction of canopy layer
% flux.albsoib   ! Direct beam albedo of ground (soil)
% flux.albsoid   ! Diffuse albedo of ground (soil)
%
% Output
% flux.swleaf    ! Leaf absorbed solar radiation (W/m2 leaf)
% flux.swveg     ! Absorbed solar radiation, vegetation (W/m2)
% flux.swvegsun  ! Absorbed solar radiation, sunlit canopy (W/m2)
% flux.swvegsha  ! Absorbed solar radiation, shaded canopy (W/m2)
% flux.swsoi     ! Absorbed solar radiation, ground (W/m2)
% flux.albcan    ! Albedo above canopy
% -----------------------------------------------------------------------

% --- Set up tridiagonal matrix

for ib = 1:params.numrad   % Process each waveband
   for p = 1:params.npts   % Process each grid point

      iv = canopy.nsoi(p);
      swup(iv) = 0;
      swdn(iv) = 0;
      for iv = canopy.nbot(p):canopy.ntop(p)
         swup(iv) = 0;
         swdn(iv) = 0;
      end

      % There are two equations for each canopy layer and the soil. The first
      % equation is the upward flux and the second equation is the downward flux. 

      m = 0; % Initialize equation index for tridiagonal matrix

      % Soil: upward flux

      iv = canopy.nsoi(p);
      m = m + 1;
      a(m) = 0;
      b(m) = 1;
      c(m) = -flux.albsoid(p,ib);
      d(m) = atmos.swskyb(p,ib) * tbcum(p,iv) * flux.albsoib(p,ib);

      % Soil: downward flux

      refld = (1 - td(p,iv+1)) * rho(p,ib);
      trand = (1 - td(p,iv+1)) * tau(p,ib) + td(p,iv+1);
      aiv = refld - trand * trand / refld;
      biv = trand / refld;

      m = m + 1;
      a(m) = -aiv;
      b(m) = 1;
      c(m) = -biv;
      d(m) = atmos.swskyb(p,ib) * tbcum(p,iv+1) * (1 - tb(p,iv+1)) * (tau(p,ib) - rho(p,ib) * biv);

      % Leaf layers, excluding top layer

      for iv = canopy.nbot(p):canopy.ntop(p)-1

         % Upward flux

         refld = (1 - td(p,iv)) * rho(p,ib);
         trand = (1 - td(p,iv)) * tau(p,ib) + td(p,iv);
         fiv = refld - trand * trand / refld;
         eiv = trand / refld;

         m = m + 1;
         a(m) = -eiv;
         b(m) = 1;
         c(m) = -fiv;
         d(m) = atmos.swskyb(p,ib) * tbcum(p,iv) * (1 - tb(p,iv)) * (rho(p,ib) - tau(p,ib) * eiv);

         % Downward flux

         refld = (1 - td(p,iv+1)) * rho(p,ib);
         trand = (1 - td(p,iv+1)) * tau(p,ib) + td(p,iv+1);
         aiv = refld - trand * trand / refld;
         biv = trand / refld;

         m = m + 1;
         a(m) = -aiv;
         b(m) = 1;
         c(m) = -biv;
         d(m) = atmos.swskyb(p,ib) * tbcum(p,iv+1) * (1 - tb(p,iv+1)) * (tau(p,ib) - rho(p,ib) * biv);

      end

      % Top canopy layer: upward flux

      iv = canopy.ntop(p);
      refld = (1 - td(p,iv)) * rho(p,ib);
      trand = (1 - td(p,iv)) * tau(p,ib) + td(p,iv);
      fiv = refld - trand * trand / refld;
      eiv = trand / refld;

      m = m + 1;
      a(m) = -eiv;
      b(m) = 1;
      c(m) = -fiv;
      d(m) = atmos.swskyb(p,ib) * tbcum(p,iv) * (1 - tb(p,iv)) * (rho(p,ib) - tau(p,ib) * eiv);

      % Top canopy layer: downward flux

      m = m + 1;
      a(m) = 0;
      b(m) = 1;
      c(m) = 0;
      d(m) = atmos.swskyd(p,ib);

      % --- Solve tridiagonal equations for fluxes

      [u] = tridiagonal_solver (a, b, c, d, m);

      % Now copy the solution (u) to the upward (swup) and downward (swdn) fluxes for each layer
      % swup - Upward diffuse solar flux above layer
      % swdn - Downward diffuse solar flux onto layer

      m = 0;

      % Soil fluxes

      iv = canopy.nsoi(p);
      m = m + 1;
      swup(iv) = u(m);
      m = m + 1;
      swdn(iv) = u(m);

      % Leaf layer fluxes

      for iv = canopy.nbot(p):canopy.ntop(p)
         m = m + 1;
         swup(iv) = u(m);
         m = m + 1;
         swdn(iv) = u(m);
      end

      % --- Compute flux densities

      % Absorbed direct beam and diffuse for ground (soil)

      iv = canopy.nsoi(p);
      direct = atmos.swskyb(p,ib) * tbcum(p,iv) * (1 - flux.albsoib(p,ib));
      diffuse = swdn(iv) * (1 - flux.albsoid(p,ib));
      flux.swsoi(p,ib) = direct + diffuse;

      % Absorbed direct beam and diffuse for each leaf layer and sum
      % for all leaf layers

      flux.swveg(p,ib) = 0;
      flux.swvegsun(p,ib) = 0;
      flux.swvegsha(p,ib) = 0;

      for iv = canopy.nbot(p):canopy.ntop(p)

         % Per unit ground area (W/m2 ground)

         direct = atmos.swskyb(p,ib) * tbcum(p,iv) * (1 - tb(p,iv)) * (1 - omega(p,ib));
         diffuse = (swdn(iv) + swup(iv-1)) * (1 - td(p,iv)) * (1 - omega(p,ib));

         % Absorbed solar radiation for shaded and sunlit portions of leaf layer
         % per unit ground area (W/m2 ground)

         sun = diffuse * flux.fracsun(p,iv) + direct;
         shade = diffuse * flux.fracsha(p,iv);

         % Convert to per unit sunlit and shaded leaf area (W/m2 leaf)

         flux.swleaf(p,iv,params.sun,ib) = sun / (flux.fracsun(p,iv) * canopy.dlai(p,iv));
         flux.swleaf(p,iv,params.sha,ib) = shade / (flux.fracsha(p,iv) * canopy.dlai(p,iv));

         % Sum fluxes over all leaf layers

         flux.swveg(p,ib) = flux.swveg(p,ib) + (direct + diffuse);
         flux.swvegsun(p,ib) = flux.swvegsun(p,ib) + sun;
         flux.swvegsha(p,ib) = flux.swvegsha(p,ib) + shade;

      end

      % --- Albedo

      incoming = atmos.swskyb(p,ib) + atmos.swskyd(p,ib);
      reflected = swup(canopy.ntop(p));
      if (incoming > 0)
         flux.albcan(p,ib) = reflected / incoming;
      else
         flux.albcan(p,ib) = 0;
      end

      % --- Conservation check

      % Total radiation balance: absorbed = incoming - outgoing

      suminc = atmos.swskyb(p,ib) + atmos.swskyd(p,ib);
      sumref = flux.albcan(p,ib) * (atmos.swskyb(p,ib) + atmos.swskyd(p,ib));
      sumabs = suminc - sumref;

      err = sumabs - (flux.swveg(p,ib) + flux.swsoi(p,ib));
      if (abs(err) > 1e-03)
         fprintf('err = %15.5f\n',err)
         fprintf('sumabs = %15.5f\n',sumabs)
         fprintf('swveg = %15.5f\n',flux.swveg(p,ib))
         fprintf('swsoi = %15.5f\n',flux.swsoi(p,ib))
         error ('NormanRadiation: Total solar conservation error')
      end

      % Sunlit and shaded absorption

      err = (flux.swvegsun(p,ib) + flux.swvegsha(p,ib)) - flux.swveg(p,ib);
      if (abs(err) > 1e-03)
         fprintf('err = %15.5f\n',err)
         fprintf('swveg = %15.5f\n',flux.swveg(p,ib))
         fprintf('swvegsun = %15.5f\n',flux.swvegsun(p,ib))
         fprintf('swvegsha = %15.5f\n',flux.swvegsha(p,ib))
         error ('NormanRadiation: Sunlit/shade solar conservation error')
      end

   end   % End grid point loop
end      % End waveband loop
