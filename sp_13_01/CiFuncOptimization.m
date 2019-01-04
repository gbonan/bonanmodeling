function [flux] = CiFuncOptimization (atmos, leaf, flux)

% Calculate leaf photosynthesis for a specified stomatal conductance.
% Then calculate Ci from the diffusion equation. 
%
% This routine uses a quadratic equation to solve for net photosynthesis (An).
% A general equation for C3 photosynthesis is:
%
%      a*(Ci - Cp)
% An = ----------- - Rd
%       e*Ci + d
%
% where:
%
% An = Net leaf photosynthesis (umol CO2/m2/s)
% Rd = Leaf respiration (umol CO2/m2/s)
% Ci = Intercellular CO2 concentration (umol/mol)
% Cp = CO2 compensation point (umol/mol)
% 
% Rubisco-limited photosynthesis (Ac)
% a  = Vcmax
% e  = 1
% d  = Kc*(1 + Oi/Ko)
%
% RuBP regeneration-limited photosynthesis (Aj)
% a = J
% e = 4
% d = 8*Cp
%
% where:
%
% Vcmax = Maximum carboxylation rate (umol/m2/s)
% Kc    = Michaelis-Menten constant for CO2 (umol/mol)
% Ko    = Michaelis-Menten constant for O2 (mmol/mol)
% Oi    = Intercellular O2 concentration (mmol/mol)
% J     = Electron transport rate (umol/m2/s)
%
% Ci is calculated from the diffusion equation:
%
%                   1.4   1.6
% An = (Ca - Ci) / (--- + ---)
%                   gb    gs
%
%            1.4   1.6
% Ci = Ca - (--- + ---)*An
%            gb    gs
%
% where:
% 
% Ca  = Atmospheric CO2 concentration (umol/mol)
% gb  = Leaf boundary layer conductance (mol H2O/m2/s)
% gs  = Leaf stomatal conductance (mol H2O/m2/s)
% 1.4 = Corrects gb for the diffusivity of CO2 compared with H2O
% 1.6 = Corrects gs for the diffusivity of CO2 compared with H2O
%
% The resulting quadratic equation is: a*An^2 + b*An + c = 0, which
% is solved for An. Correct solution is the smaller of the two roots.
%
% A similar approach is used for C4 photosynthesis.

% ------------------------------------------------------
% Input
%   atmos.o2air         ! Atmospheric O2 (mmol/mol)
%   atmos.co2air        ! Atmospheric CO2 (umol/mol)
%   leaf.c3psn          ! Photosynthetic pathway: 1 = C3. 0 = C4 plant
%   leaf.colim          ! Photosynthesis co-limitation: 0 = no. 1 = yes
%   leaf.colim_c3       ! Empirical curvature parameter for C3 co-limitation
%   leaf.colim_c4a      ! Empirical curvature parameter for C4 co-limitation
%   leaf.colim_c4b      ! Empirical curvature parameter for C4 co-limitation
%   leaf.qe_c4          ! C4: Quantum yield (mol CO2 / mol photons)
%   flux.vcmax          ! Maximum carboxylation rate (umol/m2/s)
%   flux.cp             ! CO2 compensation point (umol/mol)
%   flux.kc             ! Michaelis-Menten constant for CO2 (umol/mol)
%   flux.ko             ! Michaelis-Menten constant for O2 (mmol/mol)
%   flux.je             ! Electron transport rate (umol/m2/s)
%   flux.kp_c4          ! C4: Initial slope of CO2 response curve (mol/m2/s)
%   flux.gs             ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
%   flux.gbc            ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
%   flux.apar           ! Leaf absorbed PAR (umol photon/m2 leaf/s)
%   flux.rd             ! Leaf respiration rate (umol CO2/m2 leaf/s)
%
% Output
%   flux.ac             ! Leaf Rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.aj             ! Leaf RuBP regeneration-limited gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.ap             ! Leaf product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.ag             ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.an             ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
%   flux.cs             ! Leaf surface CO2 (umol/mol)
%   flux.ci             ! Leaf intercellular CO2 (umol/mol)
% ------------------------------------------------------

% --- Leaf conductance
% gbc has units mol CO2/m2/s
% gs has units mol H2O/m2/s
% gleaf has units mol CO2/m2/s

gleaf = 1 / (1 / flux.gbc + 1.6 / flux.gs);

% --- Gross assimilation rates

if (leaf.c3psn == 1)

   % C3: Rubisco-limited photosynthesis

   a0 = flux.vcmax;
   e0 = 1;
   d0 = flux.kc * (1 + atmos.o2air / flux.ko);

   aquad = e0 / gleaf;
   bquad = -(e0 * atmos.co2air + d0) - (a0 - e0 * flux.rd) / gleaf;
   cquad = a0 * (atmos.co2air - flux.cp) - flux.rd * (e0 * atmos.co2air + d0);

   pcoeff = [aquad bquad cquad];
   proots = roots(pcoeff);
   flux.ac = min(proots(1), proots(2)) + flux.rd;
 
   % C3: RuBP regeneration-limited photosynthesis

   a0 = flux.je;
   e0 = 4;
   d0 = 8 * flux.cp;

   aquad = e0 / gleaf;
   bquad = -(e0 * atmos.co2air + d0) - (a0 - e0 * flux.rd) / gleaf;
   cquad = a0 * (atmos.co2air - flux.cp) - flux.rd * (e0 * atmos.co2air + d0);

   pcoeff = [aquad bquad cquad];
   proots = roots(pcoeff);
   flux.aj = min(proots(1), proots(2)) + flux.rd;

   % C3: Product-limited photosynthesis

   flux.ap = 0;

else

   % C4: Rubisco-limited photosynthesis

   flux.ac = flux.vcmax;
 
   % C4: RuBP-limited photosynthesis

   flux.aj = leaf.qe_c4 * flux.apar;

   % C4: PEP carboxylase-limited (CO2-limited)

   flux.ap = flux.kp_c4 * (atmos.co2air * gleaf + flux.rd) / (gleaf + flux.kp_c4);

end

% --- Net assimilation as the minimum or co-limited rate

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

% --- Intercelluar CO2

flux.ci = atmos.co2air - flux.an / gleaf;
