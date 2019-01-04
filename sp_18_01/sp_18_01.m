% Supplemental program 18.1

% This program simulates litter decomposition and soil organic matter dynamics using DAYCENT
% as described in Figure 18.8. The program follows carbon and nitrogen pools and fluxes
% associated with the addition of a single pulse of leaf litter on the surface or root litter
% belowground. This is the litter bag protocal in the simulations of Bonan et al.
% (2013) Global Change Biology 19:957-974. This is simplified code, to illustrate the basics
% of carbon and nitrogen dynamics. The full DAYCENT model is more complex and includes
% additional processes.

% --------------
% Initialization
% --------------

% --- number of carbon pools

npool = 12;

% pools:
%  1 = metabolic litter (SRFC)
%  2 = metabolic litter (SOIL)
%  3 = structural litter (SRFC)
%  4 = structural litter (SOIL)
%  5 = coarse woody debris: fine branch
%  6 = coarse woody debris: large wood
%  7 = coarse woody debris: coarse root
%  8 = active SOM1 (SRFC)
%  9 = active SOM1 (SOIL)
% 10 = slow SOM2 (SRFC)
% 11 = slow SOM2 (SOIL)
% 12 = passive SOM3

% --- number of litter fluxes

nlitflx = 5;

%  1 = leaf
%  2 = fine root
%  3 = coarse woody debris (fine branch)
%  4 = coarse woody debris (large wood)
%  5 = coarse woody debris (coarse root)

% --- specify biome type (for mixing of surface SOM2 to soil SOM2)

grass = 1;
forest = 2;
biome = grass;

% --- site conditions

sand = 50;        % percent sand
clay = 20;        % percent clay
cdi = 0.8;        % soil temperature and moisture factor (0-1); used as a constant in this example
pH = 7.0;         % pH
O2 = 1.0;         % effect of soil anaerobic conditions on decomposition (0-1)
soilN = 3;        % soil mineral nitrogen (gN/m2); determines C/N of SOM pools

% --- effect of cultivation on decomposition (1:SOM1, 2:SOM2, 3:SOM3, 4:structural)
% These are not used here. DAYCENT has factors specific to various cultivation types.

cultfac(1) = 1;
cultfac(2) = 1;
cultfac(3) = 1;
cultfac(4) = 1;

% --- coarse woody debris chemistry

cwdlig(1) = 0.25;         % lignin fraction for fine branch
cwdlig(2) = 0.25;         % lignin fraction for large wood
cwdlig(3) = 0.25;         % lignin fraction for coarse root

% --- lignin

rsplig = 0.3;             % fraction of carbon lost as respiration (lignin)

% --- leaf and fine root litter chemistry

% generic values (need to be specified for leaf and fine root)

leaf_flig = 0.159;       % leaf litter lignin fraction
leaf_cn = 61.8;          % leaf litter C:N (gC/gN)
froot_flig = 0.349;      % fine root litter lignin fraction
froot_cn = 61.5;         % fine root C:N (gC/gN)

% chose specific litter type for litterbag study

  litter = 'TRAEf';      % leaf: Triticum aestivum
% litter = 'PIREf';      % leaf: Pinus resinosa
% litter = 'THPLf';      % leaf: Western red cedar
% litter = 'ACSAf';      % leaf: sugar maple
% litter = 'QUPRf';      % leaf: chestnut oak
% litter = 'DRGLf';      % leaf: tropical broadleaf
% litter = 'ANGEr';      % root: big bluestem grass
% litter = 'PIELr';      % root: slash pine
% litter = 'DRGLr';      % root: tropical broadleaf

% specific litter types

switch litter

   case 'TRAEf'          % leaf: Triticum aestivum
   leaf_flig = 0.162;    % leaf litter lignin fraction
   leaf_cn = 133.3;      % leaf litter C:N (gC/gN)

   case 'PIREf'          % leaf: Pinus resinosa
   leaf_flig = 0.192;    % leaf litter lignin fraction
   leaf_cn = 92.7;       % leaf litter C:N (gC/gN)

   case 'THPLf'          % leaf: Western red cedar
   leaf_flig = 0.267;    % leaf litter lignin fraction
   leaf_cn = 83.1;       % leaf litter C:N (gC/gN)

   case 'ACSAf'          % leaf: sugar maple
   leaf_flig = 0.159;    % leaf litter lignin fraction
   leaf_cn = 61.8;       % leaf litter C:N (gC/gN)

   case 'QUPRf'          % leaf: chestnut oak
   leaf_flig = 0.235;    % leaf litter lignin fraction
   leaf_cn = 50.5;       % leaf litter C:N (gC/gN)

   case 'DRGLf'          % leaf: tropical broadleaf
   leaf_flig = 0.109;    % leaf litter lignin fraction
   leaf_cn = 24.2;       % leaf litter C:N (gC/gN)

   case 'ANGEr'          % root: big bluestem grass
   froot_flig = 0.105;   % fine root litter lignin fraction
   froot_cn = 59.4;      % fine root C:N (gC/gN)

   case 'PIELr'          % root: slash pine
   froot_flig = 0.349;   % fine root litter lignin fraction
   froot_cn = 61.5;      % fine root C:N (gC/gN)

   case 'DRGLr'          % root: tropical broadleaf
   froot_flig = 0.161;   % fine root litter lignin fraction
   froot_cn = 64.6;      % fine root C:N (gC/gN)

end

% --- litter flux - U: 1=leaf, 2=fine root, 3=cwd (fine branch), 4=cwd (large wood), 5=cwd (coarse root)

% DAYCENT: fraction of large wood CWD flux that is fine branch

fbranch = 0.10;

% fluxes: gC/m2/year -> gC/m2/sec

U(1,1) = 0 / (365 * 86400);
U(2,1) = 0 / (365 * 86400);
U(3,1) = 0 / (365 * 86400) * fbranch;
U(4,1) = 0 / (365 * 86400) * (1 - fbranch);
U(5,1) = 0 / (365 * 86400);

% C/N ratios of litter fluxes

U_cn(1) = leaf_cn;
U_cn(2) = froot_cn;
U_cn(3) = 500;
U_cn(4) = 500;
U_cn(5) = 500;

% --- initial pool size
% C(i,1) = carbon for pool i (g C/m2)
% N(i,1) = nitrogen for pool i (g N/m2)

C = zeros(npool,1);
N = zeros(npool,1);

% ------------------
% Time stepping loop
% ------------------

% --- length of time step (seconds)

dt = 1800;    % 30-minutes

% --- number of timsteps per day

ntim = 86400 / dt;

% --- number of years to simulate 

nyears = 11;

% --- days per month

ndays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% --- initialize step counter

nstep = 0;

% --- counter for daily output

ncount = 0;

% --- advance time

for year = 1:nyears
   for month = 1:12
      for day = 1:ndays(month)
         for itim = 1:ntim

            nstep = nstep + 1;              % step counter
            cumday = nstep / ntim;          % cumulative day
            cumyear = nstep / (ntim * 365); % cumulative year

            if (month == 1 && day == 1 && itim == 1)
               fprintf('year = %8.1f\n',cumyear)
            end
            if (itim == 1 && day == 1)
               fprintf('day = %6.0f\n',cumday)
            end

            % =============================
            % add litter on first time step
            % =============================

            if (litter == 'TRAEf' | litter == 'PIREf' | litter == 'THPLf' | ...
                litter == 'ACSAf' | litter == 'QUPRf' | litter == 'DRGLf')
                add_leaf = 100;
                add_root = 0;
            end

            if (litter == 'ANGEr' | litter == 'PIELr' | litter == 'DRGLr')
                add_leaf = 0;
                add_root = 100;
            end

            if (nstep == 1)
               U(1,1) = add_leaf / dt;
               U(2,1) = add_root / dt;
               U(3,1) = 0;
               U(4,1) = 0;
               U(5,1) = 0;
               litterC = add_leaf + add_root;                         % save litter C input
               litterN = add_leaf / leaf_cn + add_root / froot_cn;    % save litter N input
            else
               U(1,1) = 0;
               U(2,1) = 0;
               U(3,1) = 0;
               U(4,1) = 0;
               U(5,1) = 0;
            end

            % ===============================================
            % calculate carbon fluxes and update carbon pools
            % ===============================================

            % matrix B to partition litter fluxes to each carbon pool
            % and also lignin fraction of structural litter strlig

            [strlig, B] = litter_partition_matrix (leaf_flig, leaf_cn, froot_flig, froot_cn);

            % base decomposition rate K and mixing Kmix

            [K, Kmix] = decomp_rate_base (biome, grass, forest);

            % environmental scalar xi adjusts base decomposition rate

            [K_s21, xi] = decomp_rate_scalar (cdi, pH, O2, sand, strlig, cwdlig, cultfac, K, Kmix);

            % carbon transfer matrix A
            % pathf = fractional carbon flow from pool j to pool i
            % respf = fractional respiration loss for carbon flow from pool j to pool i

            [A, pathf, respf] = carbon_transfer_matrix (npool, sand, clay, O2, strlig, cwdlig, rsplig, Kmix, K_s21);

            % calculate pool increment dC for each pool i - this
            % is the generalized calculation

            for i = 1:npool
               dC(i,1) = 0;

               % litter flux input

               for j = 1:nlitflx
                  dC(i,1) = dC(i,1) + B(i,j) * U(j,1);
               end

               % carbon transfer from pool j to pool i

               for j = 1:npool
                  if (j ~= i)
                     dC(i,1) = dC(i,1) + (1 - respf(i,j)) * pathf(i,j) * xi(j,j) * K(j,j) * C(j,1);
                  end
               end

               % carbon loss from pool i

               dC(i,1) = dC(i,1) - xi(i,i) * K(i,i) * C(i,1);

            end

            % ... or use matrix algebra: dC = B * U + A * xi * K * C

%           dC = B * U + A * xi * K * C;

            % heterotrophic respiration

            RH = 0;
            for i = 1:npool
               for j = 1:npool
                  if (j ~= i)
                     RH = RH + respf(i,j) * pathf(i,j) * xi(j,j) * K(j,j) * C(j,1);
                  end
               end
            end

            % =====================================
            % calculate N fluxes and update N pools
            % =====================================

            % specify C/N ratio of SOM pools based on soil mineral N

            cnrat = @(x,cnmin,cnmax,ncrit) max(min(((1-x/ncrit)*(cnmax-cnmin)+cnmin),cnmax),cnmin);

            for i = 1:7
               som_cn(i) = 0;
            end
            som_cn( 8) = cnrat(soilN, 10, 20, 1);
            som_cn( 9) = cnrat(soilN,  8, 18, 2);
            som_cn(10) = cnrat(soilN, 12, 15, 2);
            som_cn(11) = cnrat(soilN, 12, 40, 2);
            som_cn(12) = cnrat(soilN,  6, 20, 2);

            % N input from litter flux j to pool i. note the
            % special case for leaf and fine root litter
            % to the structural pool (which has C/N = 200).

            N_in_total = U(1,1) / U_cn(1);
            N_to_struc = B(3,1) * U(1,1) / 200;
            N_to_metab = N_in_total - N_to_struc;
            N_litter_flux(1) = N_to_metab;
            N_litter_flux(3) = N_to_struc;

            N_in_total = U(2,1) / U_cn(2);
            N_to_struc = B(4,2) * U(2,1) / 200;
            N_to_metab = N_in_total - N_to_struc;
            N_litter_flux(2) = N_to_metab;
            N_litter_flux(4) = N_to_struc;

            N_litter_flux(5) = B(5,3) * U(3,1) / U_cn(3);
            N_litter_flux(6) = B(6,4) * U(4,1) / U_cn(4);
            N_litter_flux(7) = B(7,5) * U(5,1) / U_cn(5);

            N_litter_flux(8) = 0;
            N_litter_flux(9) = 0;
            N_litter_flux(10) = 0;
            N_litter_flux(11) = 0;
            N_litter_flux(12) = 0;

            % N flow in carbon transfer among pools

            for i = 1:npool
               N_transfer(i) = 0;

               % N transfer from donor pool j to receiver pool i

               for j = 1:npool
                  if (j ~= i)
                     N_transfer(i) = N_transfer(i) + pathf(i,j) * xi(j,j) * K(j,j) * N(j,1);
                  end
               end

               % N loss from pool i

               N_transfer(i) = N_transfer(i) - xi(i,i) * K(i,i) * N(i,1);

            end

            % mineral N flux in carbon transfer from donor pool j to SOM receiver pool i

            for i = 1:12
               N_mineral_flux(i) = 0;
            end

            for i = 8:12
               for j = 1:npool

                  % carbon flow from pool j to pool i

                  tcflow = pathf(i,j) * xi(j,j) * K(j,j) * C(j,1);

                  if (tcflow > 0)

                     % C/N ratio of donor pool

                     cn_donor = C(j,1) / N(j,1);

                     % C/N ratio of receiver SOM pool

                     cn_receiver = som_cn(i);

                     % special case for SOM2(SRFC) -> SOM2(SOIL); there is no immobilization
                     % or mineralization because this is a mixing flux

                     if (j == 10 & i == 11)
                        cn_receiver = cn_donor;
                     end

                     % Note the special case when lignin decomposes to SOM1 and SOM2. In the
                     % DAYCENT code used by Bonan et al. (2013), the receiver C/N ratios are:
                     % SOM1: cn_receiver = cnrat(soilN, 10, 20, 1)
                     % SOM2: cn_receiver = cnrat(soilN, 12, 15, 2)
                     % and do not vary between SRFC and SOIL pools. This means that the target
                     % C/N for SOIL pools is the same as for SRFC pools.

                     % special case when lignin decomposes to SOM1(SOIL): use values for SOM1(SRFC)

                     if ((j == 4 & i == 9) | (j == 7 & i == 9))
                        cn_receiver = cnrat(soilN, 10, 20, 1);
                     end

                     % special case when lignin decomposes to SOM2(SOIL): use values for SOM2(SRFC)

                     if ((j == 4 & i == 11) | (j == 7 & i == 11))
                        cn_receiver = cnrat(soilN, 12, 15, 2);
                     end

                     % mineral N flux for flow from donor to receiver (+ = immobilization)

                     if (j ~= i)
                        N_flux = tcflow * ((1 - respf(i,j)) / cn_receiver - 1/cn_donor);
                     end

                     % sum N flux to pool i

                     N_mineral_flux(i) = N_mineral_flux(i) + N_flux;

                  end
               end
            end

            % sum N fluxes for each pool i

            for i = 1:npool
               dN(i,1) = N_litter_flux(i) + N_transfer(i) + N_mineral_flux(i);
            end

            % ============
            % update pools
            % ============

            % carbon pools

            for i = 1:npool
               C(i,1) = C(i,1) + dC(i,1) * dt;
            end

            litC = C(1,1) + C(2,1) + C(3,1) + C(4,1);             % litter: metabolic + structural
            cwdC = C(5,1) + C(6,1) + C(7,1);                      % coarse woody debris
            somC = C(8,1) + C(9,1) + C(10,1) + C(11,1) + C(12,1); % soil organic matter
            totC = litC + cwdC + somC;                            % total carbon

            % balance check

            dCtot = 0;
            for i = 1:npool
               dCtot = dCtot + dC(i,1);
            end

            Utot = 0;
            for i = 1:nlitflx
               Utot = Utot + U(i,1);
            end

            err = (dCtot - (Utot - RH)) * dt;
            if (abs(err) > 1e-12)
               fprintf('err = %15.5f\n',err)
               error ('CARBON BALANCE ERROR')
            end

            % nitrogen pools

            totN_old = 0;
            for i = 1:npool
               totN_old = totN_old + N(i,1);
            end

            for i = 1:npool
               N(i,1) = N(i,1) + dN(i,1) * dt;
            end

            litN = N(1,1) + N(2,1) + N(3,1) + N(4,1);             % litter: metabolic + structural
            cwdN = N(5,1) + N(6,1) + N(7,1);                      % coarse woody debris
            somN = N(8,1) + N(9,1) + N(10,1) + N(11,1) + N(12,1); % soil organic matter
            totN = litN + cwdN + somN;                            % total nitrogen

            % balance check

            Utot = 0;
            Nmin_tot = 0;
            for i = 1:npool
               Utot = Utot + N_litter_flux(i);
               Nmin_tot = Nmin_tot + N_mineral_flux(i);
            end

            dNtot = totN - totN_old;
            err = dNtot - (Utot + Nmin_tot) * dt;
            if (abs(err) > 1e-12)
               fprintf('err = %15.5f\n',err)
               error ('NITROGEN BALANCE ERROR')
            end

            % save daily output for graphing

            if (itim == 1)
               ncount = ncount + 1;
               x1(ncount) = cumyear;
               y1(ncount) = totC;
               y2(ncount) = totN;
           end

         end
      end
   end
end

% ----------------------
% graph data
% ----------------------

% plot C and N vs time

plot(x1,y1/litterC,'g-',x1,y2/litterN,'r-')
xlabel('Year')
ylabel('Fraction of initial mass')
legend('carbon','nitrogen','Location','best')
title(litter)

% plot N vs C

% plot(y1/litterC,y2/litterN)
% set(gca, 'XDir','reverse')
% xlabel('Mass remaining (fraction)')
% ylabel('Fraction of initial nitrogen')
% title(litter)

% ----------------------------------
% Write formated output to text file
% ----------------------------------

data = [x1; y1; y2];
fileID = fopen('data.txt','w');
fprintf(fileID,'%15s %15s %15s\n','year','carbon','nitrogen');
fprintf(fileID,'%15.5f %15.5f %15.5f\n',data);
