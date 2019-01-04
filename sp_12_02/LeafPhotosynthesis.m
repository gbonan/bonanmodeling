function [flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux)

% Calculate leaf photosynthesis using one of two methods:
% (1) Calculate leaf photosynthesis and stomatal conductance
% by solving for the value of Ci that satisfies the
% metabolic, stomatal constraint, and diffusion equations.
% This is used with the Ball-Berry style stomatal models.
% (2) Calculate leaf photosynthesis for a specified stomatal
% conductance. Then calculate Ci from the diffusion equation.
% This is used with the WUE stomatal optimization.

% ------------------------------------------------------
% Input
%   physcon.tfrz        ! Freezing point of water (K)
%   physcon.rgas        ! Universal gas constant (J/K/mol)
%   atmos.co2air        ! Atmospheric CO2 (umol/mol)
%   atmos.eair          ! Vapor pressure of air (Pa)
%   leaf.gstyp          ! Stomatal conductance: 0 = Medlyn. 1 = Ball-Berry. 2 = WUE optimization
%   leaf.c3psn          ! Photosynthetic pathway: 1 = C3. 0 = C4 plant
%   leaf.vcmax25        ! Maximum carboxylation rate at 25C (umol/m2/s)
%   leaf.jmax25         ! Maximum electron transport rate at 25C (umol/m2/s)
%   leaf.rd25           ! Leaf respiration rate at 25C (umol CO2/m2/s)
%   leaf.kc25           ! Michaelis-Menten constant for CO2 at 25C (umol/mol)
%   leaf.ko25           ! Michaelis-Menten constant for O2 at 25C (mmol/mol)
%   leaf.cp25           ! CO2 compensation point at 25C (umol/mol)
%   leaf.kcha           ! Activation energy for Kc (J/mol)
%   leaf.koha           ! Activation energy for Ko (J/mol)
%   leaf.cpha           ! Activation energy for Cp (J/mol)
%   leaf.vcmaxha        ! Activation energy for Vcmax (J/mol)
%   leaf.jmaxha         ! Activation energy for Jmax (J/mol)
%   leaf.rdha           ! Activation energy for Rd (J/mol)
%   leaf.vcmaxhd        ! Deactivation energy for Vcmax (J/mol)
%   leaf.jmaxhd         ! Deactivation energy for Jmax (J/mol)
%   leaf.rdhd           ! Deactivation energy for Rd (J/mol)
%   leaf.vcmaxse        ! Entropy term for Vcmax (J/mol/K)
%   leaf.jmaxse         ! Entropy term for Jmax (J/mol/K)
%   leaf.rdse           ! Entropy term for Rd (J/mol/K)
%   leaf.vcmaxc         ! Vcmax scaling factor for high temperature inhibition (25 C = 1.0)
%   leaf.jmaxc          ! Jmax scaling factor for high temperature inhibition (25 C = 1.0)
%   leaf.rdc            ! Rd scaling factor for high temperature inhibition (25 C = 1.0)
%   leaf.phi_psii       ! Quantum yield of PS II
%   leaf.theta_j        ! Empirical curvature parameter for electron transport rate
%   leaf.kp25_c4        ! C4: Initial slope of CO2 response curve at 25C (mol/m2/s)
%   leaf.g0             ! Ball-Berry minimum leaf conductance (mol H2O/m2/s)
%   leaf.g1             ! Ball-Berry slope of conductance-photosynthesis relationship
%   flux.gbv            ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
%   flux.gbc            ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
%   flux.apar           ! Leaf absorbed PAR (umol photon/m2 leaf/s)
%   flux.tleaf          ! Leaf temperature (K)
%
% Input or output (depending on method)
%   flux.gs             ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
%
% Output
%   flux.vcmax          ! Maximum carboxylation rate (umol/m2/s)
%   flux.jmax           ! Maximum electron transport rate (umol/m2/s)
%   flux.cp             ! CO2 compensation point (umol/mol)
%   flux.kc             ! Michaelis-Menten constant for CO2 (umol/mol)
%   flux.ko             ! Michaelis-Menten constant for O2 (mmol/mol)
%   flux.je             ! Electron transport rate (umol/m2/s)
%   flux.kp_c4          ! C4: Initial slope of CO2 response curve (mol/m2/s)
%   flux.ac             ! Leaf Rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.aj             ! Leaf RuBP regeneration-limited gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.ap             ! Leaf product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.ag             ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.an             ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
%   flux.rd             ! Leaf respiration rate (umol CO2/m2 leaf/s)
%   flux.cs             ! Leaf surface CO2 (umol/mol)
%   flux.ci             ! Leaf intercellular CO2 (umol/mol)
%   flux.hs             ! Leaf fractional humidity at surface (dimensionless)
%   flux.vpd            ! Leaf vapor pressure deficit at surface (Pa)
% ------------------------------------------------------

% --- Adjust photosynthetic parameters for temperature

if (leaf.c3psn == 1)

   % C3 temperature response

   ft = @(tl, ha) exp(ha/(physcon.rgas*(physcon.tfrz+25)) * (1-(physcon.tfrz+25)/tl));
   fth = @(tl, hd, se, fc) fc / (1 + exp((-hd+se*tl)/(physcon.rgas*tl)));

   flux.kc = leaf.kc25 * ft(flux.tleaf, leaf.kcha);
   flux.ko = leaf.ko25 * ft(flux.tleaf, leaf.koha);
   flux.cp = leaf.cp25 * ft(flux.tleaf, leaf.cpha);

   t1 = ft(flux.tleaf, leaf.vcmaxha);
   t2 = fth(flux.tleaf, leaf.vcmaxhd, leaf.vcmaxse, leaf.vcmaxc);
   flux.vcmax = leaf.vcmax25 * t1 * t2;

   t1 = ft(flux.tleaf, leaf.jmaxha);
   t2 = fth(flux.tleaf, leaf.jmaxhd, leaf.jmaxse, leaf.jmaxc);
   flux.jmax = leaf.jmax25 * t1 * t2;

   t1 = ft(flux.tleaf, leaf.rdha);
   t2 = fth(flux.tleaf, leaf.rdhd, leaf.rdse, leaf.rdc);
   flux.rd = leaf.rd25 * t1 * t2;

   flux.kp_c4 = 0;

elseif (leaf.c3psn == 0)

   % C4 temperature response

   t1 = 2^((flux.tleaf-(physcon.tfrz+25)) / 10);
   t2 = 1 + exp(0.2*((physcon.tfrz+15)-flux.tleaf));
   t3 = 1 + exp(0.3*(flux.tleaf-(physcon.tfrz+40)));
   flux.vcmax = leaf.vcmax25 * t1 / (t2 * t3);

   t3 = 1 + exp(1.3*(flux.tleaf-(physcon.tfrz+55)));
   flux.rd = leaf.rd25  * t1 / t3;

   flux.kp_c4 = leaf.kp25_c4 * t1;

   flux.kc = 0;
   flux.ko = 0;
   flux.cp = 0;
   flux.jmax = 0;
   flux.je = 0;

end

% --- Electron transport rate for C3 plants

% Solve the polynomial: aquad*Je^2 + bquad*Je + cquad = 0
% for Je. Correct solution is the smallest of the two roots.

if (leaf.c3psn == 1)
   qabs = 0.5 * leaf.phi_psii * flux.apar;
   aquad = leaf.theta_j;
   bquad = -(qabs + flux.jmax);
   cquad = qabs * flux.jmax;
   pcoeff = [aquad bquad cquad];
   proots = roots(pcoeff);
   flux.je = min(proots(1), proots(2));
end

% --- Ci calculation

if (leaf.gstyp <= 1)

   % Initial estimates for Ci

   if (leaf.c3psn == 1)
      ci0 = 0.7 * atmos.co2air;
   elseif (leaf.c3psn == 0)
      ci0 = 0.4 * atmos.co2air;
   end
   ci1 = ci0 * 0.99;

   % Solve for Ci: Use CiFunc to iterate photosynthesis calculations
   % until the change in Ci is < tol. Ci has units umol/mol

   tol = 0.1;                 % Accuracy tolerance for Ci (umol/mol)
   func_name = 'CiFunc';      % The function name

   [flux, dummy] = hybrid_root (func_name, physcon, atmos, leaf, flux, ci0, ci1, tol);
   flux.ci = dummy;

elseif (leaf.gstyp == 2)

   % Calculate photosynthesis for a specified stomatal conductance

   [flux] = CiFuncOptimization (atmos, leaf, flux);

end

% --- Relative humidity and vapor pressure at leaf surface

[esat, desat] = satvap ((flux.tleaf-physcon.tfrz));
flux.hs = (flux.gbv * atmos.eair + flux.gs * esat) / ((flux.gbv + flux.gs) * esat);
flux.vpd = max(esat - flux.hs*esat, 0.1);

% --- Make sure iterative solution is correct

if (flux.gs < 0)
   error ('LeafPhotosynthesis: negative stomatal conductance')
end

% Compare with Ball-Berry model. The solution blows up with low eair. In input
% data, eair should be > 0.05*esat to ensure that hs does not go to zero.

if (leaf.gstyp == 1)
   gs_err = leaf.g1 * max(flux.an, 0) * flux.hs / flux.cs + leaf.g0;
   if (abs(flux.gs-gs_err)*1e06 > 1e-04)
      fprintf('gs = %15.4f\n', flux.gs)
      fprintf('gs_err = %15.4f\n', gs_err)
      error ('LeafPhotosynthesis: failed Ball-Berry error check')
   end
end

% Compare with Medlyn model. The solutions blows up with vpd = 0. The
% quadratic calcuation of gsw in CiFunc constrains vpd > 50 Pa, so this
% comparison is only valid for those conditions.

if (leaf.gstyp == 0)
   if ((esat - atmos.eair) > 50)
      gs_err = 1.6 * (1 + leaf.g1 / sqrt(flux.vpd*0.001)) * max(flux.an, 0) / flux.cs + leaf.g0;
      if (abs(flux.gs-gs_err)*1e06 > 1e-04)
         fprintf('gs = %15.4f\n', flux.gs)
         fprintf('gs_err = %15.4f\n', gs_err)
         error ('LeafPhotosynthesis: failed Medlyn error check')
      end
   end
end

% Compare with diffusion equation: An = (ca - ci) * gleaf

an_err = (atmos.co2air - flux.ci) / (1 / flux.gbc + 1.6 / flux.gs);
if (flux.an > 0 & abs(flux.an-an_err) > 0.01)
   fprintf('An = %15.4f\n', flux.an)
   fprintf('An_err = %15.4f\n', an_err)
   error ('LeafPhotosynthesis: failed diffusion error check')
end
