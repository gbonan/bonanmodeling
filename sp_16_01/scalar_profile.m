function [fluxvar] = scalar_profile (dt, p, physcon, forcvar, surfvar, leafvar, soilvar, fluxvar)

% Compute scalar profiles for temperature and water vapor using an implicit
% solution. The boundary condition is the above-canopy scalar values at the
% reference height. The vegetation and ground temperature and fluxes are
% calculated as part of the implicit solution.

% -----------------------------------------------------------------------------------
% Input
%   dt                  ! Model time step (s)
%   p                   ! Index for grid point to process
%   physcon.tfrz        ! Freezing point of water (K)
%   physcon.hvap        ! Latent heat of evaporation (J/kg)
%   physcon.mmh2o       ! Molecular mass of water (kg/mol)
%   forcvar.rhomol      ! Molar density (mol/m3)
%   forcvar.cpair       ! Specific heat of air at constant pressure (J/mol/K)
%   forcvar.pref        ! Atmospheric pressure (Pa)
%   forcvar.thref       ! Potential temperature at reference height (K)
%   forcvar.qref        ! Water vapor at reference height (mol/mol)
%   surfvar.nlev        ! Index for top level
%   surfvar.ntop        ! Index for top leaf layer
%   surfvar.nsoi        ! First canopy layer is soil
%   surfvar.zw          ! Canopy height at layer interfaces (m)
%   surfvar.dpai        ! Layer plant area index (m2/m2)
%   surfvar.fwet        ! Fraction of plant area index that is wet
%   surfvar.fdry        ! Fraction of plant area index that is green and dry
%   surfvar.fracsun     ! Sunlit fraction of canopy layer
%   surfvar.fracsha     ! Shaded fraction of canopy layer
%   leafvar.nleaf       ! Number of leaf types (sunlit and shaded)
%   leafvar.isun        ! Sunlit leaf index
%   leafvar.isha        ! Shaded leaf index
%   leafvar.gbh         ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
%   leafvar.gbv         ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
%   leafvar.gs          ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
%   leafvar.cpleaf      ! Leaf heat capacity (J/m2 leaf/K)
%   leafvar.rnleaf      ! Leaf net radiation (W/m2 leaf)
%   soilvar.tk          ! Soil thermal conductivity (W/m/K)
%   soilvar.dz          ! Soil layer depth (m)
%   soilvar.tsoi        ! Soil temperature (K)
%   soilvar.resis       ! Soil evaporative resistance (s/m)
%   soilvar.rhg         ! Relative humidity of airspace at soil surface (fraction)
%   fluxvar.rnsoi       ! Net radiation at ground (W/m2)
%   fluxvar.tair_old    ! Air temperature profile for previous timestep (K)
%   fluxvar.qair_old    ! Water vapor profile for previous timestep (mol/mol)
%   fluxvar.tveg_old    ! Vegetation temperature profile for previous timestep (K)
%   fluxvar.ga_prof     ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
% Output
%   fluxvar.tair        ! Air temperature profile (K)
%   fluxvar.qair        ! Water vapor profile (mol/mol)
%   fluxvar.tveg        ! Vegetation temperature profile (K)
%   fluxvar.shsoi       ! Ground sensible heat flux, ground (W/m2)
%   fluxvar.etsoi       ! Ground evaporation flux (mol H2O/m2/s)
%   fluxvar.gsoi        ! Soil heat flux (W/m2)
%   fluxvar.tg          ! Soil surface temperature (K)
%   fluxvar.shveg       ! Leaf sensible heat flux (W/m2 leaf)
%   fluxvar.etveg       ! Leaf evapotranspiration flux (mol H2O/m2 leaf/s)
%   fluxvar.stveg       ! Leaf storage heat flux (W/m2 leaf)
%   fluxvar.shair       ! Canopy air sensible heat flux (W/m2)
%   fluxvar.etair       ! Canopy air water vapor flux (mol H2O/m2/s)
%   fluxvar.stair       ! Canopy air storage heat flux (W/m2)
% -----------------------------------------------------------------------------------

% Latent heat of vaporization (J/mol)

lambda = physcon.hvap * physcon.mmh2o;

% --------------------------------------------------------------------------
% Terms for ground temperature, which is calculated from the energy balance:
% Rn0 - H0 - lambda*E0 - G = 0
% and is rewritten as:
% T0 = alpha0*T(1) + beta0*q(1) + delta0
% --------------------------------------------------------------------------

% array index for ground

ic = surfvar.nsoi(p);

% Saturation vapor pressure and derivative at ground temperature

[esat, desat] = satvap (fluxvar.tair_old(p,ic)-physcon.tfrz); % Pa
qsat0 = esat / forcvar.pref(p);                               % Pa -> mol/mol
dqsat0 = desat / forcvar.pref(p);                             % Pa -> mol/mol

% Total soil-to-air conductance for water vapor (mol H2O/m2/s)

gsw = (1 / soilvar.resis(p)) * forcvar.rhomol(p);             % soil: s/m -> mol H2O/m2/s
gsa = fluxvar.ga_prof(p,ic);                                  % aerodynamic
gs0 = gsa * gsw / (gsa + gsw);                                % total conductance

% Terms for soil heat flux

c02 = soilvar.tk(p) / soilvar.dz(p);
c01 = -c02 * soilvar.tsoi(p);

% Coefficients for ground temperature

den = forcvar.cpair(p) * fluxvar.ga_prof(p,ic) + lambda * soilvar.rhg(p) * gs0 * dqsat0 + c02;
alpha0 = forcvar.cpair(p) * fluxvar.ga_prof(p,ic) / den;
beta0 = lambda * gs0 / den;
delta0 = (fluxvar.rnsoi(p) - lambda * soilvar.rhg(p) * gs0 * (qsat0 - dqsat0 * fluxvar.tair_old(p,ic)) - c01) / den;

% ---------------------------------------------------------------------
% alpha, beta, delta coefficients for leaf temperature:
%
% Tlsun(i) = alpha_sun(i)*T(i) + beta_sun(i)*q(i) + delta_sun(i)
% Tlsha(i) = alpha_sha(i)*T(i) + beta_sha(i)*q(i) + delta_sha(i)
%
% and leaf terms needed for scalar conservations equations. These
% are defined for all layers but the leaf terms only exist for
% the layers with leaves.
% ---------------------------------------------------------------------

for ic = surfvar.nsoi(p)+1:surfvar.nlev(p)

   if (surfvar.dpai(p,ic) > 0)

      % Calculate terms for sunlit and shaded leaves

      for il = 1:leafvar.nleaf

         % Leaf conductances (here these are per unit leaf area)

         gleaf_sh(ic,il) = 2 * leafvar.gbh(p,ic,il);
         gleaf = leafvar.gs(p,ic,il) * leafvar.gbv(p,ic,il) / (leafvar.gs(p,ic,il) + leafvar.gbv(p,ic,il));
         gleaf_et(ic,il) = gleaf * surfvar.fdry(p,ic) + leafvar.gbv(p,ic,il) * surfvar.fwet(p,ic);

         % Heat capacity of leaves

         heatcap(ic,il) = leafvar.cpleaf(p,ic);

         % Available energy: net radiation

         avail_energy(ic,il) = leafvar.rnleaf(p,ic,il);

         % Saturation vapor pressure and derivative for leaf temperature at time n: Pa -> mol/mol

         [esat, desat] = satvap (fluxvar.tveg_old(p,ic,il)-physcon.tfrz);
         qsat(ic,il) = esat / forcvar.pref(p);
         dqsat(ic,il) = desat / forcvar.pref(p);

         % Term for linearized vapor pressure at leaf temperature:
         % qsat(tveg) = qsat(tveg_old) + dqsat * (tveg - tveg_old)
         % Here qsat_term contains the terms with tveg_old

         qsat_term(ic,il) = qsat(ic,il) - dqsat(ic,il) * fluxvar.tveg_old(p,ic,il);

         % alpha, beta, delta coefficients for leaf temperature

         den = heatcap(ic,il) / dt + gleaf_sh(ic,il) * forcvar.cpair(p) + gleaf_et(ic,il) * lambda * dqsat(ic,il);
         alpha(ic,il) = gleaf_sh(ic,il) * forcvar.cpair(p) / den;
         beta(ic,il) = gleaf_et(ic,il) * lambda / den;
         delta(ic,il) = avail_energy(ic,il) / den ...
                      - lambda * gleaf_et(ic,il) * qsat_term(ic,il) / den ...
                      + heatcap(ic,il) / dt * fluxvar.tveg_old(p,ic,il) / den;

         % Now scale flux terms for plant area

         if (il == leafvar.isun)
            pai = surfvar.fracsun(p,ic) * surfvar.dpai(p,ic);
         elseif (il == leafvar.isha)
            pai = surfvar.fracsha(p,ic) * surfvar.dpai(p,ic);
         end

         gleaf_sh(ic,il) = gleaf_sh(ic,il) * pai;
         gleaf_et(ic,il) = gleaf_et(ic,il) * pai;
         heatcap(ic,il) = heatcap(ic,il) * pai;
         avail_energy(ic,il) = avail_energy(ic,il) * pai;

      end

   else

      % Zero out terms

      for il = 1:leafvar.nleaf
         gleaf_sh(ic,il) = 0;
         gleaf_et(ic,il) = 0;
         heatcap(ic,il) = 0;
         avail_energy(ic,il) = 0;
         qsat(ic,il) = 0;
         dqsat(ic,il) = 0;
         qsat_term(ic,il) = 0;
         alpha(ic,il) = 0;
         beta(ic,il) = 0;
         delta(ic,il) = 0;
      end

   end
end

% ---------------------------------------------------------------------
% a,b,c,d coefficients for air temperature:
% a1(i)*T(i-1) + b11(i)*T(i) + b12(i)*q(i) + c1(i)*T(i+1) = d1(i)
%
% a,b,c,d coefficients for water vapor (mole fraction):
% a2(i)*q(i-1) + b21(i)*T(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)
% ---------------------------------------------------------------------

for ic = surfvar.nsoi(p)+1:surfvar.nlev(p)

   % Storage term

   rho_dz_over_dt = forcvar.rhomol(p) * (surfvar.zw(p,ic) - surfvar.zw(p,ic-1)) / dt;

   % a1,b11,b12,c1,d1 coefficients for air temperature

   a1(ic) = -fluxvar.ga_prof(p,ic-1);
   b11(ic) = rho_dz_over_dt + fluxvar.ga_prof(p,ic-1) + fluxvar.ga_prof(p,ic) ...
           + gleaf_sh(ic,leafvar.isun) * (1 - alpha(ic,leafvar.isun)) ...
           + gleaf_sh(ic,leafvar.isha) * (1 - alpha(ic,leafvar.isha));
   b12(ic) = -gleaf_sh(ic,leafvar.isun) * beta(ic,leafvar.isun) - gleaf_sh(ic,leafvar.isha) * beta(ic,leafvar.isha);
   c1(ic) = -fluxvar.ga_prof(p,ic);
   d1(ic) = rho_dz_over_dt * fluxvar.tair_old(p,ic) + gleaf_sh(ic,leafvar.isun) * delta(ic,leafvar.isun) ...
          + gleaf_sh(ic,leafvar.isha) * delta(ic,leafvar.isha);

   % Special case for top layer

   if (ic == surfvar.nlev(p))
      c1(ic) = 0;
      d1(ic) = d1(ic) + fluxvar.ga_prof(p,ic) * forcvar.thref(p);
   end

   % Special case for first canopy layer (i.e., immediately above the ground)

   if (ic == surfvar.nsoi(p)+1)
      a1(ic) = 0;
      b11(ic) = b11(ic) - fluxvar.ga_prof(p,surfvar.nsoi(p)) * alpha0;
      b12(ic) = b12(ic) - fluxvar.ga_prof(p,surfvar.nsoi(p)) * beta0;
      d1(ic) = d1(ic) + fluxvar.ga_prof(p,surfvar.nsoi(p)) * delta0;
   end

   % a2,b21,b22,c2,d2 coefficients for water vapor (mole fraction)

   if (ic == surfvar.nsoi(p)+1)
      ga_prof_ic_minus_one = gs0;
   else
      ga_prof_ic_minus_one = fluxvar.ga_prof(p,ic-1);
   end

   a2(ic) = -ga_prof_ic_minus_one;
   b21(ic) = -gleaf_et(ic,leafvar.isun) * dqsat(ic,leafvar.isun) * alpha(ic,leafvar.isun) ...
             -gleaf_et(ic,leafvar.isha) * dqsat(ic,leafvar.isha) * alpha(ic,leafvar.isha);
   b22(ic) = rho_dz_over_dt + ga_prof_ic_minus_one + fluxvar.ga_prof(p,ic) ...
           + gleaf_et(ic,leafvar.isun) * (1 - dqsat(ic,leafvar.isun) * beta(ic,leafvar.isun)) ...
           + gleaf_et(ic,leafvar.isha) * (1 - dqsat(ic,leafvar.isha) * beta(ic,leafvar.isha));
   c2(ic) = -fluxvar.ga_prof(p,ic);
   d2(ic) = rho_dz_over_dt * fluxvar.qair_old(p,ic) ...
          + gleaf_et(ic,leafvar.isun) * (dqsat(ic,leafvar.isun) * delta(ic,leafvar.isun) + qsat_term(ic,leafvar.isun)) ...
          + gleaf_et(ic,leafvar.isha) * (dqsat(ic,leafvar.isha) * delta(ic,leafvar.isha) + qsat_term(ic,leafvar.isha));

   % Special case for top layer

   if (ic == surfvar.nlev(p))
      c2(ic) = 0;
      d2(ic) = d2(ic) + fluxvar.ga_prof(p,ic) * forcvar.qref(p);
   end

   % Special case for first canopy layer (i.e., immediately above the ground)

   if (ic == surfvar.nsoi(p)+1)
      a2(ic) = 0;
      b21(ic) = b21(ic) - gs0 * soilvar.rhg(p) * dqsat0 * alpha0;
      b22(ic) = b22(ic) - gs0 * soilvar.rhg(p) * dqsat0 * beta0;
      d2(ic) = d2(ic) + gs0 * soilvar.rhg(p) * (qsat0 + dqsat0 * (delta0 - fluxvar.tair_old(p,surfvar.nsoi(p))));
   end

end

% ---------------------------------------------------------------------
% Solve for air temperature and water vapor (mole fraction):
%
% a1(i)*T(i-1) + b11(i)*T(i) + b12(i)*q(i) + c1(i)*T(i+1) = d1(i)
% a2(i)*q(i-1) + b21(i)*T(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)
%
% The solution rewrites these equations so that:
% T(i) = f1(i) - e11(i)*T(i+1) - e12(i)*q(i+1) 
% q(i) = f2(i) - e21(i)*T(i+1) - e22(i)*q(i+1) 
% ---------------------------------------------------------------------

ic = surfvar.nsoi(p);
e11(ic) = 0;
e12(ic) = 0;
e21(ic) = 0;
e22(ic) = 0;
f1(ic) = 0;
f2(ic) = 0;

for ic = surfvar.nsoi(p)+1:surfvar.nlev(p)

   % The matrix to invert is:
   %
   % B(i)- A(i)*E(i-1)
   %
   % which is a 2x2 matrix. The
   % elements in the 2x2 matrix are: 
   %
   %                     | a b |
   % B(i)- A(i)*E(i-1) = | c d |
   %
   % Calculate these elements (denoted by ainv,binv,
   % cinv,dinv) and the determinant of the matrix.

   ainv = b11(ic) - a1(ic) * e11(ic-1);
   binv = b12(ic) - a1(ic) * e12(ic-1);
   cinv = b21(ic) - a2(ic) * e21(ic-1);
   dinv = b22(ic) - a2(ic) * e22(ic-1);
   det = ainv * dinv - binv * cinv;

   % E(i) = [B(i) - A(i)*E(i-1)]^(-1) * C(i)

   e11(ic) = dinv * c1(ic) / det;
   e12(ic) = -binv * c2(ic) / det;
   e21(ic) = -cinv * c1(ic) / det;
   e22(ic) = ainv * c2(ic) / det;

   % F(i) = [B(i) - A(i)*E(i-1)]^(-1) * [D(i) - A(i)*F(i-1)]

   f1(ic) =  (dinv*(d1(ic) - a1(ic)*f1(ic-1)) - binv*(d2(ic) - a2(ic)*f2(ic-1))) / det;
   f2(ic) = (-cinv*(d1(ic) - a1(ic)*f1(ic-1)) + ainv*(d2(ic) - a2(ic)*f2(ic-1))) / det;

end

% Top layer

ic = surfvar.nlev(p);
fluxvar.tair(p,ic) = f1(ic);
fluxvar.qair(p,ic) = f2(ic);

% Layers through to bottom of canopy

for ic = surfvar.nlev(p)-1: -1: surfvar.nsoi(p)+1
   fluxvar.tair(p,ic) = f1(ic) - e11(ic)*fluxvar.tair(p,ic+1) - e12(ic)*fluxvar.qair(p,ic+1);
   fluxvar.qair(p,ic) = f2(ic) - e21(ic)*fluxvar.tair(p,ic+1) - e22(ic)*fluxvar.qair(p,ic+1);
end

% Ground

ic = surfvar.nsoi(p);
fluxvar.tg(p) = alpha0 * fluxvar.tair(p,ic+1) + beta0 * fluxvar.qair(p,ic+1) + delta0;
fluxvar.tair(p,ic) = fluxvar.tg(p);
fluxvar.qair(p,ic) = soilvar.rhg(p) * (qsat0 + dqsat0 * (fluxvar.tair(p,ic) - fluxvar.tair_old(p,ic)));

% ---------------------------------------------------------------------
% Calculate leaf temperature:
% Tlsun(i) = alpha_sun(i)*T(i) + beta_sun(i)*q(i) + delta_sun(i)
% Tlsha(i) = alpha_sha(i)*T(i) + beta_sha(i)*q(i) + delta_sha(i)
% ---------------------------------------------------------------------

for ic = surfvar.nsoi(p)+1:surfvar.nlev(p)
   fluxvar.tveg(p,ic,leafvar.isun) = alpha(ic,leafvar.isun)*fluxvar.tair(p,ic) ...
                                   + beta(ic,leafvar.isun)*fluxvar.qair(p,ic) + delta(ic,leafvar.isun);
   fluxvar.tveg(p,ic,leafvar.isha) = alpha(ic,leafvar.isha)*fluxvar.tair(p,ic) ...
                                   + beta(ic,leafvar.isha)*fluxvar.qair(p,ic) + delta(ic,leafvar.isha);

   % Special checks for no plant area in layer

   if (surfvar.dpai(p,ic) > 0)
      pai = surfvar.fracsun(p,ic) * surfvar.dpai(p,ic);
      if (pai == 0)
         fluxvar.tveg(p,ic,leafvar.isun) = fluxvar.tair(p,ic);
      end
      pai = surfvar.fracsha(p,ic) * surfvar.dpai(p,ic);
      if (pai == 0)
         fluxvar.tveg(p,ic,leafvar.isha) = fluxvar.tair(p,ic);
      end
   else
      fluxvar.tveg(p,ic,leafvar.isun) = 0;
      fluxvar.tveg(p,ic,leafvar.isha) = 0;
   end
end

ic = surfvar.nsoi(p);
fluxvar.tveg(p,ic,leafvar.isun) = 0;
fluxvar.tveg(p,ic,leafvar.isha) = 0;

% --------------------
% Ground source fluxes
% --------------------

ic = surfvar.nsoi(p);
fluxvar.shsoi(p) = -forcvar.cpair(p) * (fluxvar.tair(p,ic+1) - fluxvar.tair(p,ic)) * fluxvar.ga_prof(p,ic);
fluxvar.etsoi(p) = -(fluxvar.qair(p,ic+1) - fluxvar.qair(p,ic)) * gs0;
fluxvar.gsoi(p) = c01 + c02 * fluxvar.tg(p);

% ------------------------
% Vegetation source fluxes
% ------------------------

for ic = surfvar.nsoi(p)+1:surfvar.nlev(p)
   fluxvar.shveg(p,ic) = 0;
   fluxvar.etveg(p,ic) = 0;
   fluxvar.stveg(p,ic) = 0;
   for il = 1:leafvar.nleaf
      fluxvar.shveg(p,ic) = fluxvar.shveg(p,ic) ...
                          + forcvar.cpair(p) * (fluxvar.tveg(p,ic,il) - fluxvar.tair(p,ic)) * gleaf_sh(ic,il);
      fluxvar.etveg(p,ic) = fluxvar.etveg(p,ic) ...
                          + (qsat(ic,il) + dqsat(ic,il) * (fluxvar.tveg(p,ic,il) - fluxvar.tveg_old(p,ic,il)) - fluxvar.qair(p,ic)) ...
                          * gleaf_et(ic,il);
      fluxvar.stveg(p,ic) = fluxvar.stveg(p,ic) ...
                          + heatcap(ic,il) * (fluxvar.tveg(p,ic,il) - fluxvar.tveg_old(p,ic,il)) / dt;
   end
end

% ------------------------------------------------------------
% Vertical sensible heat and water vapor fluxes between layers
% ------------------------------------------------------------

for ic = surfvar.nsoi(p)+1:surfvar.nlev(p)-1
   fluxvar.shair(p,ic) = -forcvar.cpair(p) * (fluxvar.tair(p,ic+1) - fluxvar.tair(p,ic)) * fluxvar.ga_prof(p,ic);
   fluxvar.etair(p,ic) = -(fluxvar.qair(p,ic+1) - fluxvar.qair(p,ic)) * fluxvar.ga_prof(p,ic);
end

ic = surfvar.nlev(p);
fluxvar.shair(p,ic) = -forcvar.cpair(p) * (forcvar.thref(p) - fluxvar.tair(p,ic)) * fluxvar.ga_prof(p,ic);
fluxvar.etair(p,ic) = -(forcvar.qref(p) - fluxvar.qair(p,ic)) * fluxvar.ga_prof(p,ic);

% --------------------------------------------------------------------------
% Canopy air storage flux (W/m2) and its sensible heat and water vapor terms
% --------------------------------------------------------------------------

for ic = surfvar.nsoi(p)+1:surfvar.nlev(p)
   dz_over_dt = (surfvar.zw(p,ic) - surfvar.zw(p,ic-1)) / dt;
   storage_sh(ic) = forcvar.rhomol(p) * forcvar.cpair(p) * (fluxvar.tair(p,ic) - fluxvar.tair_old(p,ic)) * dz_over_dt;
   storage_et(ic) = forcvar.rhomol(p) * (fluxvar.qair(p,ic) - fluxvar.qair_old(p,ic)) * dz_over_dt;
   fluxvar.stair(p,ic) = storage_sh(ic) + storage_et(ic) * lambda;
end

% -------------------
% Conservation checks
% -------------------

% Ground source fluxes energy balance

err = fluxvar.rnsoi(p) - fluxvar.shsoi(p) - lambda * fluxvar.etsoi(p) - fluxvar.gsoi(p);
if (abs(err) > 0.001)
   error ('ScalarProfile: Ground temperature energy balance error')
end

% Vegetation source fluxes energy balance

for ic = surfvar.nsoi(p)+1:surfvar.nlev(p)
   err = avail_energy(ic,leafvar.isun) + avail_energy(ic,leafvar.isha) ...
       - fluxvar.shveg(p,ic) - lambda * fluxvar.etveg(p,ic) - fluxvar.stveg(p,ic);
   if (abs(err) > 0.001)
      error ('ScalarProfile: Leaf temperature energy balance error')
   end
end

% Flux conservation at each layer. Note special case for first canopy layer.

for ic = surfvar.nsoi(p)+1:surfvar.nlev(p)

   if (ic == surfvar.nsoi(p)+1)
      err = storage_sh(ic) - (fluxvar.shsoi(p) + fluxvar.shveg(p,ic) - fluxvar.shair(p,ic));
   else
      err = storage_sh(ic) - (fluxvar.shair(p,ic-1) + fluxvar.shveg(p,ic) - fluxvar.shair(p,ic));
   end
   if (abs(err) > 0.001)
      error ('ScalarProfile: Sensible heat layer conservation error')
   end

   if (ic == surfvar.nsoi(p)+1)
      err = storage_et(ic) - (fluxvar.etsoi(p) + fluxvar.etveg(p,ic) - fluxvar.etair(p,ic));
   else
      err = storage_et(ic) - (fluxvar.etair(p,ic-1) + fluxvar.etveg(p,ic) - fluxvar.etair(p,ic));
   end
   err = err * lambda;
   if (abs(err) > 0.001)
      err ('ScalarProfile: Latent heat layer conservation error')
   end

end

% Flux conservation for canopy sensible heat. This is to check canopy
% conservation equation (so the sum is to ntop not nlev).

sum_src = fluxvar.shsoi(p);
sum_storage = 0;
for ic = surfvar.nsoi(p)+1:surfvar.ntop(p)
   sum_src = sum_src + fluxvar.shveg(p,ic);
   sum_storage = sum_storage + storage_sh(ic);
end

err = (sum_src - sum_storage) - fluxvar.shair(p,surfvar.ntop(p));
if (abs(err) > 0.001)
   error ('ScalarProfile: Sensible heat canopy conservation error')
end

% Flux conservation for canopy latent heat. This is to check canopy
% conservation equation (so the sum is to ntop not nlev).

sum_src = fluxvar.etsoi(p);
sum_storage = 0;
for ic = surfvar.nsoi(p)+1:surfvar.ntop(p)
   sum_src = sum_src + fluxvar.etveg(p,ic);
   sum_storage = sum_storage + storage_et(ic);
end

err = (sum_src - sum_storage) - fluxvar.etair(p,surfvar.ntop(p));
err = err * lambda;
if (abs(err) > 0.001)
   error ('ScalarProfile: Latent heat canopy conservation error')
end
