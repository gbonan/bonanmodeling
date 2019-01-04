function [A, pathf, respf] = carbon_transfer_matrix (npool, sand, clay, O2, strlig, cwdlig, rsplig, Kmix, K_s21)

% Carbon transfer matrix

% ---------------------------------------------------------------------------------------------
% Input
%   npool       ! number of carbon pools
%   sand        ! percent sand
%   clay        ! percent clay
%   O2          ! effect of soil anaerobic conditions on decomposition (0-1)
%   strlig      ! lignin fraction: (1) SRFC and (2) SOIL structural litter (g lignin/g biomass)
%   cwdlig      ! lignin fraction: (1) fine branch; (2) large wood; (3) coarse root
%   rsplig      ! fraction of carbon lost as respiration (lignin)
%   Kmix        ! base mixing rate: SOM2(SRFC) -> SOM2(SOIL), 1/sec
%   K_s21       ! rate constant: total loss from SOM2(SRFC), 1/sec
%
% Output
%   A           ! carbon transfer matrix
%   pathf       ! fractional carbon flow from pool j to pool i
%   respf       ! fractional respiration loss for carbon flow from pool j to pool i
% ---------------------------------------------------------------------------------------------

% --- zero out arrays

for i = 1:npool
   for j = 1:npool
      pathf(i,j) = 0;
      respf(i,j) = 0;
      A(i,j) = 0;
   end
end

% --- anaerobic factor

fanerb = 1 + 5 * (1 - O2);

% --- pathf(i,j) = fractional carbon flow from pool j to pool i

pathf(8,1) = 1;
pathf(9,2) = 1;
pathf(8,3) = 1 - strlig(1);
pathf(10,3) = strlig(1);
pathf(9,4) = 1 - strlig(2);
pathf(11,4) = strlig(2);

pathf(8,5) = 1 - cwdlig(1);
pathf(10,5) = cwdlig(1);
pathf(8,6) = 1 - cwdlig(2);
pathf(10,6) = cwdlig(2);
pathf(9,7) = 1 - cwdlig(3);
pathf(11,7) = cwdlig(3);

pathf(10,8) = 1;
pathf(12,9) = (0.003 + 0.032 * clay/100) * fanerb;
pathf(11,9) = 1 - pathf(12,9);
pathf(11,10) = Kmix / K_s21;
pathf(8,10) = 1 - pathf(11,10);
pathf(12,11) = (0.003 + 0.009 * clay/100) * fanerb;
pathf(9,11) = 1 - pathf(12,11);
pathf(9,12) = 1;

% --- respf(i,j) = fractional respiration loss for carbon flow from pool j to pool i

respf(8,1) = 0.55;
respf(9,2) = 0.55;
respf(8,3) = 0.45;
respf(10,3) = rsplig;
respf(9,4) = 0.55;
respf(11,4) = rsplig;

respf(8,5) = 0.45;
respf(10,5) = rsplig;
respf(8,6) = 0.45;
respf(10,6) = rsplig;
respf(9,7) = 0.55;
respf(11,7) = rsplig;

respf(10,8) = 0.60;
respf(12,9) = 0;
respf(11,9) = (0.17 + 0.68 * sand/100) / pathf(11,9);
respf(11,10) = 0;
respf(8,10) = 0.55;
respf(12,11) = 0;
respf(9,11) = 0.55 / pathf(9,11);
respf(9,12) = 0.55 * O2;

% --- carbon transfer matrix: A(i,j) = fractional carbon flow from pool j that enters pool i

for i = 1:npool
   A(i,i) = -1;
end

for i = 1:npool
   for j = 1:npool
      if (j ~= i)
         A(i,j) = pathf(i,j) * (1 - respf(i,j));
      end
   end
end
