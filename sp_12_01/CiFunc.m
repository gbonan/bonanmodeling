function [flux, ci_dif] = CiFunc (physcon, atmos, leaf, flux, ci_val)

% Calculate leaf photosynthesis and stomatal conductance for a specified Ci
% (ci_val). Then calculate a new Ci from the diffusion equation. This function
% returns a value ci_dif = 0 when Ci has converged to the value that satisfies
% the metabolic, stomatal constraint, and diffusion equations.

% ------------------------------------------------------
% Input
%   physcon.tfrz        ! Freezing point of water (K)
%   atmos.o2air         ! Atmospheric O2 (mmol/mol)
%   atmos.co2air        ! Atmospheric CO2 (umol/mol)
%   atmos.eair          ! Vapor pressure of air (Pa)
%   leaf.c3psn          ! Photosynthetic pathway: 1 = C3. 0 = C4 plant
%   leaf.gstyp          ! Stomatal conductance: 0 = Medlyn. 1 = Ball-Berry. 2 = WUE optimization
%   leaf.colim          ! Photosynthesis co-limitation: 0 = no. 1 = yes
%   leaf.colim_c3       ! Empirical curvature parameter for C3 co-limitation
%   leaf.colim_c4a      ! Empirical curvature parameter for C4 co-limitation
%   leaf.colim_c4b      ! Empirical curvature parameter for C4 co-limitation
%   leaf.qe_c4          ! C4: Quantum yield (mol CO2 / mol photons)
%   leaf.g0             ! Ball-Berry minimum leaf conductance (mol H2O/m2/s)
%   leaf.g1             ! Ball-Berry slope of conductance-photosynthesis relationship
%   flux.vcmax          ! Maximum carboxylation rate (umol/m2/s)
%   flux.cp             ! CO2 compensation point (umol/mol)
%   flux.kc             ! Michaelis-Menten constant for CO2 (umol/mol)
%   flux.ko             ! Michaelis-Menten constant for O2 (mmol/mol)
%   flux.je             ! Electron transport rate (umol/m2/s)
%   flux.kp_c4          ! C4: Initial slope of CO2 response curve (mol/m2/s)
%   flux.rd             ! Leaf respiration rate (umol CO2/m2 leaf/s)
%   flux.gbv            ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
%   flux.gbc            ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
%   flux.apar           ! Leaf absorbed PAR (umol photon/m2 leaf/s)
%   flux.tleaf          ! Leaf temperature (K)
%   ci_val              ! Input value for Ci (umol/mol)
%
% Output
%   flux.ac             ! Leaf Rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.aj             ! Leaf RuBP regeneration-limited gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.ap             ! Leaf product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.ag             ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.an             ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
%   flux.cs             ! Leaf surface CO2 (umol/mol)
%   flux.gs             ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
%   ci_dif              ! Difference in Ci
% ------------------------------------------------------

% --- Metabolic (demand-based) photosynthetic rate

if (leaf.c3psn == 1)

   % C3: Rubisco-limited photosynthesis
   flux.ac = flux.vcmax * max(ci_val - flux.cp, 0) / (ci_val + flux.kc * (1 + atmos.o2air / flux.ko));
 
   % C3: RuBP regeneration-limited photosynthesis
   flux.aj = flux.je * max(ci_val - flux.cp, 0) / (4 * ci_val + 8 * flux.cp);

   % C3: Product-limited photosynthesis (do not use)
   flux.ap = 0;

else

   % C4: Rubisco-limited photosynthesis
   flux.ac = flux.vcmax;
 
   % C4: RuBP regeneration-limited photosynthesis
   flux.aj = leaf.qe_c4 * flux.apar;

   % C4: PEP carboxylase-limited (CO2-limited)
   flux.ap = flux.kp_c4 * max(ci_val, 0);

end

% --- Net photosynthesis as the minimum or co-limited rate

if (leaf.colim == 1) % Use co-limitation

   % First co-limit Ac and Aj. Ai is the intermediate co-limited photosynthesis
   % rate found by solving the polynomial: aquad*Ai^2 + bquad*Ai + cquad = 0 for Ai.
   % Correct solution is the smallest of the two roots.

   if (leaf.c3psn == 1)
      aquad = leaf.colim_c3;
   else
      aquad = leaf.colim_c4a;
   end
   bquad = -(flux.ac + flux.aj);
   cquad = flux.ac * flux.aj;
   pcoeff = [aquad bquad cquad];
   proots = roots(pcoeff);
   ai = min(proots(1), proots(2));

   % Now co-limit again using Ap, but only for C4 plants. Solve the polynomial:
   % aquad*Ag^2 + bquad*Ag + cquad = 0 for Ag. Correct solution is the smallest
   % of the two roots. Ignore the product-limited rate For C3 plants.

   if (leaf.c3psn == 0)
      aquad = leaf.colim_c4b;
      bquad = -(ai + flux.ap);
      cquad = ai * flux.ap;
      pcoeff = [aquad bquad cquad];
      proots = roots(pcoeff);
      flux.ag = min(proots(1), proots(2));
   else
      flux.ag = ai;
   end

elseif (leaf.colim == 0) % No co-limitation

   if (leaf.c3psn == 1)
      flux.ag = min(flux.ac, flux.aj); % C3
   else
      flux.ag = min(flux.ac, flux.aj, flux.ap); % C4
   end

end

% Prevent photosynthesis from ever being negative

flux.ac = max(flux.ac, 0);
flux.aj = max(flux.aj, 0);
flux.ap = max(flux.ap, 0);
flux.ag = max(flux.ag, 0);

% Net CO2 uptake

flux.an = flux.ag - flux.rd;

% --- CO2 at leaf surface

flux.cs = atmos.co2air - flux.an / flux.gbc;
flux.cs = max(flux.cs, 1);

% --- Stomatal constraint function

% Saturation vapor pressure at leaf temperature

[esat, desat] = satvap ((flux.tleaf-physcon.tfrz));

% Ball-Berry stomatal conductance is a quadratic equation
% for gs given An: aquad*gs^2 + bquad*gs + cquad = 0. Correct
% solution is the larger of the two roots. This solution is
% valid for An >= 0. With An <= 0, gs = g0.

if (leaf.gstyp == 1)
   term = flux.an / flux.cs;
   if (flux.an > 0)
      aquad = 1;
      bquad = flux.gbv - leaf.g0 - leaf.g1 * term;
      cquad = -flux.gbv * (leaf.g0 + leaf.g1 * term * atmos.eair / esat);
      pcoeff = [aquad bquad cquad];
      proots = roots(pcoeff);
      flux.gs = max(proots(1), proots(2));
   else
      flux.gs = leaf.g0;
   end
end

% Quadratic equation for Medlyn stomatal conductance

if (leaf.gstyp == 0)
   vpd = max((esat - atmos.eair), 50) * 0.001;
   term = 1.6 * flux.an / flux.cs;
   if (flux.an > 0)
      aquad = 1;
      bquad = -(2 * (leaf.g0 + term) + (leaf.g1 * term)^2 / (flux.gbv * vpd));
      cquad = leaf.g0 * leaf.g0 + (2 * leaf.g0 + term * (1 - leaf.g1 * leaf.g1 / vpd)) * term;
      pcoeff = [aquad bquad cquad];
      proots = roots(pcoeff);
      flux.gs = max(proots(1), proots(2));
   else
      flux.gs = leaf.g0;
   end
end

% --- Diffusion (supply-based) photosynthetic rate

% Leaf CO2 conductance (mol CO2/m2/s)

gleaf = 1 / (1 / flux.gbc + 1.6 / flux.gs);

% Calculate Ci from the diffusion rate

cinew = atmos.co2air - flux.an / gleaf;

% --- Return the difference between the current Ci and the new Ci

if (flux.an >= 0)
   ci_dif = cinew - ci_val;
else
   ci_dif = 0;
end
