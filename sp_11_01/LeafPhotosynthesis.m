function [flux] = LeafPhotosynthesis (physcon, atmos, leaf, flux)

% Calculate leaf photosynthesis given a known Ci

% ------------------------------------------------------
% Input
%   physcon.tfrz        ! Freezing point of water (K)
%   physcon.rgas        ! Universal gas constant (J/K/mol)
%   atmos.co2air        ! Atmospheric CO2 (umol/mol)
%   atmos.o2air         ! Atmospheric O2 (mmol/mol)
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
%   leaf.colim          ! Photosynthesis co-limitation: 0 = no. 1 = yes
%   leaf.colim_c3       ! Empirical curvature parameter for C3 co-limitation
%   flux.apar           ! Leaf absorbed PAR (umol photon/m2 leaf/s)
%   flux.tleaf          ! Leaf temperature (K)
%
% Output
%   flux.vcmax          ! Maximum carboxylation rate (umol/m2/s)
%   flux.jmax           ! Maximum electron transport rate (umol/m2/s)
%   flux.cp             ! CO2 compensation point (umol/mol)
%   flux.kc             ! Michaelis-Menten constant for CO2 (umol/mol)
%   flux.ko             ! Michaelis-Menten constant for O2 (mmol/mol)
%   flux.je             ! Electron transport rate (umol/m2/s)
%   flux.ac             ! Leaf Rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.aj             ! Leaf RuBP regeneration-limited gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.ag             ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
%   flux.an             ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
%   flux.rd             ! Leaf respiration rate (umol CO2/m2 leaf/s)
%   flux.ci             ! Leaf intercellular CO2 (umol/mol)
% ------------------------------------------------------

% --- Adjust photosynthetic parameters for temperature

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

% --- Electron transport rate for C3 plants

% Solve the polynomial: aquad*je^2 + bquad*je + cquad = 0
% for Je. Correct solution is the smallest of the two roots.

qabs = 0.5 * leaf.phi_psii * flux.apar;
aquad = leaf.theta_j;
bquad = -(qabs + flux.jmax);
cquad = qabs * flux.jmax;
pcoeff = [aquad bquad cquad];
proots = roots(pcoeff);
flux.je = min(proots(1), proots(2));

% --- Specify Ci

flux.ci = 0.7 * atmos.co2air;

% --- C3: Rubisco-limited photosynthesis

flux.ac = flux.vcmax * max(flux.ci - flux.cp, 0) / (flux.ci + flux.kc * (1 + atmos.o2air / flux.ko));
 
% --- C3: RuBP regeneration-limited photosynthesis

flux.aj = flux.je * max(flux.ci - flux.cp, 0) / (4 * flux.ci + 8 * flux.cp);

% --- Photosynthesis is the minimum or co-limited rate

if (leaf.colim == 1) % Use co-limitation

   % Co-limit Ac and Aj. Ag is the co-limited photosynthesis rate found
   % by solving the polynomial: aquad*Ag^2 + bquad*Ag + cquad = 0 for Ag.
   % Correct solution is the smallest of the two roots.

   aquad = leaf.colim_c3;
   bquad = -(flux.ac + flux.aj);
   cquad = flux.ac * flux.aj;
   pcoeff = [aquad bquad cquad];
   proots = roots(pcoeff);
   flux.ag = min(proots(1), proots(2));

elseif (leaf.colim == 0) % No co-limitation

   flux.ag = min(flux.ac, flux.aj);

end

% Net CO2 uptake

flux.an = flux.ag - flux.rd;
