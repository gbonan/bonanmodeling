function [K, Kmix] = decomp_rate_base (biome, grass, forest)

% -----------------------------------------------------
% Base decomposition rate for each carbon pool (per sec)
% -----------------------------------------------------

% convert rate from per year to per second

K(1,1) = 8.0 / (365 * 86400);
K(2,2) = 18.5 / (365 * 86400);
K(3,3) = 2.0 / (365 * 86400);
K(4,4) = 4.9 / (365 * 86400);
K(5,5) = 1.5 / (365 * 86400);
K(6,6) = 0.02 / (365 * 86400);
K(7,7) = 0.1 / (365 * 86400);
K(8,8) = 6.0 / (365 * 86400);
K(9,9) = 11.0 / (365 * 86400);
K(10,10) = 0.08 / (365 * 86400);
K(11,11) = 0.4 / (365 * 86400);
K(12,12) = 0.0033 / (365 * 86400);

% mixing of som2c(SRFC) to som2c(SOIL)
% annual turnover time for grass/crop: 4 years
% annual turnover time for forest: 10 years

if (biome == grass)
   Kmix = 0.25 / (365 * 86400);
elseif (biome == forest)
   Kmix = 0.10 / (365 * 86400);
else
   Kmix = 0;
end
