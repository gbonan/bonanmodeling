function [leaf] = LeafPhysiologyParams (params, physcon, leaf)

% ------------------------------------------------------
% Input
%   params.vis          ! Waveband index for visible radiation
%   params.nir          ! Waveband index for near-infrared radiation
%   physcon.tfrz        ! Freezing point of water (K)
%   physcon.rgas        ! Universal gas constant (J/K/mol)
%   leaf.c3psn          ! Photosynthetic pathway: 1 = C3. 0 = C4 plant
%   leaf.gstyp          ! Stomatal conductance: 0 = Medlyn. 1 = Ball-Berry. 2 = WUE optimization
%
% Output
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
%   leaf.colim_c3       ! Empirical curvature parameter for C3 co-limitation
%   leaf.colim_c4a      ! Empirical curvature parameter for C4 co-limitation
%   leaf.colim_c4b      ! Empirical curvature parameter for C4 co-limitation
%   leaf.qe_c4          ! C4: Quantum yield (mol CO2 / mol photons)
%   leaf.kp25_c4        ! C4: Initial slope of CO2 response curve at 25C (mol/m2/s)
%   leaf.g0             ! Minimum leaf conductance (mol H2O/m2/s)
%   leaf.g1             ! Slope of conductance-photosynthesis relationship
%   leaf.dleaf          ! Leaf dimension (m)
%   leaf.emiss          ! Leaf emissivity
%   leaf.rho            ! Leaf reflectance for visible (vis) and near-infrared (nir) wavebands
%   leaf.tau            ! Leaf transmittance for visible (vis) and near-infrared (nir) wavebands
%   leaf.iota           ! Stomatal efficiency (umol CO2/ mol H2O)
% ------------------------------------------------------

% --- Vcmax and other parameters (at 25C)

if (leaf.c3psn == 1)
   leaf.vcmax25 = 60;
   leaf.jmax25 = 1.67 * leaf.vcmax25;
   leaf.kp25_c4 = 0;
   leaf.rd25 = 0.015 * leaf.vcmax25;
else
   leaf.vcmax25 = 40;
   leaf.jmax25 = 0;
   leaf.kp25_c4 = 0.02 * leaf.vcmax25;
   leaf.rd25 = 0.025 * leaf.vcmax25;
end

% --- Kc, Ko, Cp at 25C

leaf.kc25 = 404.9;
leaf.ko25 = 278.4;
leaf.cp25 = 42.75;

% --- Activation energy

leaf.kcha = 79430;
leaf.koha = 36380;
leaf.cpha = 37830;
leaf.rdha = 46390;
leaf.vcmaxha = 65330;
leaf.jmaxha  = 43540;

% --- High temperature deactivation

% Deactivation energy (J/mol)

leaf.rdhd = 150000;
leaf.vcmaxhd = 150000;
leaf.jmaxhd  = 150000;

% Entropy term (J/mol/K)

leaf.rdse = 490;
leaf.vcmaxse = 490;
leaf.jmaxse  = 490;

% Scaling factors for high temperature inhibition (25 C = 1.0).
% The factor "c" scales the deactivation to a value of 1.0 at 25C.

fth25 = @(hd, se) 1 + exp((-hd + se*(physcon.tfrz+25)) / (physcon.rgas*(physcon.tfrz+25)));
leaf.vcmaxc = fth25 (leaf.vcmaxhd, leaf.vcmaxse);
leaf.jmaxc  = fth25 (leaf.jmaxhd, leaf.jmaxse);
leaf.rdc    = fth25 (leaf.rdhd, leaf.rdse);

% --- C3 parameters

% Quantum yield of PS II

leaf.phi_psii = 0.85;

% Empirical curvature parameter for electron transport rate

leaf.theta_j = 0.90;

% Empirical curvature parameter for C3 co-limitation

leaf.colim_c3 = 0.98;

% Empirical curvature parameters for C4 co-limitation

leaf.colim_c4a = 0.80;
leaf.colim_c4b = 0.95;

% --- C4: Quantum yield (mol CO2 / mol photons)

leaf.qe_c4 = 0.05;

% --- Stomatal conductance parameters

if (leaf.c3psn == 1)

   if (leaf.gstyp == 1)
      leaf.g0 = 0.01;       % Ball-Berry minimum leaf conductance (mol H2O/m2/s)
      leaf.g1 = 9.0;        % Ball-Berry slope of conductance-photosynthesis relationship
   elseif (leaf.gstyp == 0)
      leaf.g0 = 0.0;        % Medlyn minimum leaf conductance (mol H2O/m2/s)
      leaf.g1 = 4.45;       % Medlyn slope of conductance-photosynthesis relationship
   end

else

   if (leaf.gstyp == 1)
      leaf.g0 = 0.04;       % Ball-Berry minimum leaf conductance (mol H2O/m2/s)
      leaf.g1 = 4.0;        % Ball-Berry slope of conductance-photosynthesis relationship
   elseif (leaf.gstyp == 0)
      leaf.g0 = 0.0;        % Medlyn minimum leaf conductance (mol H2O/m2/s)
      leaf.g1 = 1.62;       % Medlyn slope of conductance-photosynthesis relationship
   end

end

% --- Stomatal efficiency for optimization (An/E; umol CO2/ mol H2O)

if (leaf.gstyp == 2)
   leaf.iota = 750;
end

% --- Leaf dimension (m)

leaf.dleaf = 0.05;

% --- Leaf emissivity

leaf.emiss = 0.98;

% --- Leaf reflectance and transmittance: visible and near-infrared wavebands

leaf.rho(params.vis) = 0.10;
leaf.tau(params.vis) = 0.10;
leaf.rho(params.nir) = 0.40;
leaf.tau(params.nir) = 0.40;
