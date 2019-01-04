function [K_s21, xi] = decomp_rate_scalar (cdi, pH, O2, sand, strlig, cwdlig, cultfac, K, Kmix)

% Environmental scalar for each carbon pool adjusts base
% rate for soil temperature and soil moisture scalars
% (cdi) and additionally pH, lignin, texture, anaerobic,
% and cultivation

% ----------------------------------------------------------------
% Input
%   cdi         ! soil temperature and moisture scalar (0-1)
%   pH          ! pH
%   O2          ! effect of soil anaerobic conditions on decomposition (0-1)
%   sand        ! percent sand
%   strlig      ! lignin fraction: (1) SRFC and (2) SOIL structural litter (g lignin/g biomass)
%   cwdlig      ! lignin fraction: (1) fine branch; (2) large wood; (3) coarse root
%   cultfac     ! effect of cultivation on decomposition (1:SOM1, 2:SOM2, 3:SOM3, 4:structural)
%   K           ! base decomposition rate (1/sec)
%   Kmix        ! base mixing rate: SOM2(SRFC) -> SOM2(SOIL), 1/sec
%
% Output
%   K_s21       ! rate constant: total loss from SOM2(SRFC), 1/sec
%   xi          ! environmental scalar
% ----------------------------------------------------------------

% function for pH effect

pH_factor = @(x,a,b,c,d) b + (c / pi) * atan(d * (x - a) * pi);

% ==========================================
% metabolic litter (SRFC)
% ==========================================

pHeff = pH_factor(pH, 4.8, 0.5, 1.14, 0.7);
pHeff = max(min(pHeff, 1), 0);
xi(1,1) = cdi * pHeff;

% ==========================================
% metabolic litter (SOIL)
% ==========================================

pHeff = pH_factor(pH, 4.8, 0.5, 1.14, 0.7);
pHeff = max(min(pHeff, 1), 0);
xi(2,2) = cdi * pHeff * O2;

% ==========================================
% structural litter (SRFC)
% ==========================================

pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
pHeff = max(min(pHeff, 1), 0);
xi(3,3) = cdi * pHeff * exp(-3*strlig(1));

% ==========================================
% structural litter (SOIL)
% ==========================================

pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
pHeff = max(min(pHeff, 1), 0);
xi(4,4) = cdi * pHeff * exp(-3*strlig(2)) * O2 * cultfac(4);

% ==========================================
% coarse woody debris: fine branch
% ==========================================

pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
pHeff = max(min(pHeff, 1), 0);
xi(5,5) = cdi * pHeff * exp(-3*cwdlig(1));

% ==========================================
% coarse woody debris: large wood
% ==========================================

pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
pHeff = max(min(pHeff, 1), 0);
xi(6,6) = cdi * pHeff * exp(-3*cwdlig(2));

% ==========================================
% coarse woody debris: coarse root
% ==========================================

pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
pHeff = max(min(pHeff, 1), 0);
xi(7,7) = cdi * pHeff * exp(-3*cwdlig(3)) * O2;

% ==========================================
% active soil organic matter: SOM1 (SRFC)
% ==========================================

pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
pHeff = max(min(pHeff, 1), 0);
xi(8,8) = cdi * pHeff;

% ==========================================
% active soil organic matter: SOM1 (SOIL)
% ==========================================

pHeff = pH_factor(pH, 4.8, 0.5, 1.14, 0.7);
pHeff = max(min(pHeff, 1), 0);
texture = 0.25 + 0.75 * (sand/100);
xi(9,9) = cdi * pHeff * O2 * texture * cultfac(1);

% ==========================================
% slow soil organic matter: SOM2 (SRFC)
% ==========================================

% som2(SRFC) -> som1(SRFC)

pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
pHeff = max(min(pHeff, 1), 0);

K_s21_to_s11 = K(10,10) * pHeff;

% som2(SRFC) -> som2(SOIL): mixing

K_s21_to_s22 = Kmix;

% total loss from som2(SRFC)

K_s21 = K_s21_to_s11 + K_s21_to_s22;

% effective environmental scalar
xi(10,10) = cdi * (K_s21 / K(10,10));

% ==========================================
% slow soil organic matter: SOM2 (SOIL)
% ==========================================

pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
pHeff = max(min(pHeff, 1), 0);
xi(11,11) = cdi * pHeff * O2 * cultfac(2);

% ==========================================
% passive soil organic matter: SOM3
% ==========================================

pHeff = pH_factor(pH, 3.0, 0.5, 1.1, 0.7);
pHeff = max(min(pHeff, 1), 0);
xi(12,12) = cdi * pHeff * O2 * cultfac(3);
