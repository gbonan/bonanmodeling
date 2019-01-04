function [strlig, B] = litter_partition_matrix (leaf_flig, leaf_cn, froot_flig, froot_cn)

% Matrix to partition litter fluxes to each carbon pool

% ---------------------------------------------------------------------------------------------
% Input
%   leaf_flig   ! leaf litter lignin fraction
%   leaf_cn     ! leaf litter C:N (gC/gN)
%   froot_flig  ! fine root litter lignin fraction
%   froot_cn    ! fine root C:N (gC/gN)
%
% Output
%   strlig      ! lignin fraction: (1) SRFC and (2) SOIL structural litter (g lignin/g biomass)
%   B           ! litter flux partitioning matrix
% ---------------------------------------------------------------------------------------------

% --- leaf

% fraction of plant residue that is nitrogen (g N / g biomass)

fnit = 1 / (leaf_cn * 2.5);

% lignin fraction of plant residue (g lignin / g biomass)

flig = leaf_flig;

% lignin/nitrogen ratio of plant residue

rlig2n = flig / fnit;

% fraction of plant residue that goes to metabolic litter pool

fmet = 0.85 - 0.013 * rlig2n;

% make sure the fraction of residue which is lignin is not greater than the
% fraction which goes to structural

if (flig > (1 - fmet))
   fmet = (1 - flig);
end

% minimum metabolic fraction

if (fmet < 0.20) 
   fmet = 0.20;
end

% adjust lignin content of structural litter pool. fligst is the fraction of
% incoming structural residue that is lignin

fligst = min(flig/(1 - fmet), 1);

% DAYCENT adjusts lignin content of structural litter pool for new material
% that is added. That is not needed in this program because only one type
% of leaf litter is added. If different types of leaf litter are added,
% need to adjust structural litter lignin

strlig(1) = fligst;

% B(i,j) = litter flux partitioning matrix (litter flux j -> pool i)

B( 1,1) = fmet;
B( 2,1) = 0;
B( 3,1) = 1 - fmet;
B( 4,1) = 0;
B( 5,1) = 0;
B( 6,1) = 0;
B( 7,1) = 0;
B( 8,1) = 0;
B( 9,1) = 0;
B(10,1) = 0;
B(11,1) = 0;
B(12,1) = 0;

% --- fine root

% fraction of plant residue that is nitrogen (g N / g biomass)

fnit = 1 / (froot_cn * 2.5);

% lignin fraction of plant residue (g lignin / g biomass)

flig = froot_flig;

% lignin/nitrogen ratio of plant residue

rlig2n = flig / fnit;

% fraction of plant residue that goes to metabolic litter pool

fmet = 0.85 - 0.013 * rlig2n;

% make sure the fraction of residue which is lignin is not greater than the
% fraction which goes to structural

if (flig > (1 - fmet))
   fmet = (1 - flig);
end

% minimum metabolic fraction

if (fmet < 0.20)
   fmet = 0.20;
end

% adjust lignin content of structural litter pool. fligst is the fraction of
% incoming structural residue that is lignin

fligst = min(flig/(1 - fmet), 1);

% DAYCENT adjusts lignin content of structural litter pool for new material
% that is added. That is not needed in this program because only one type
% of fine root litter is added. If different types of fine root litter are
% added, need to adjust structural litter lignin

strlig(2) = fligst;

% B(i,j) = litter flux partitioning matrix (litter flux j -> pool i)

B( 1,2) = 0;
B( 2,2) = fmet;
B( 3,2) = 0;
B( 4,2) = 1 - fmet;
B( 5,2) = 0;
B( 6,2) = 0;
B( 7,2) = 0;
B( 8,2) = 0;
B( 9,2) = 0;
B(10,2) = 0;
B(11,2) = 0;
B(12,2) = 0;

% --- cwd (fine branch)

B( 1,3) = 0;
B( 2,3) = 0;
B( 3,3) = 0;
B( 4,3) = 0;
B( 5,3) = 1;
B( 6,3) = 0;
B( 7,3) = 0;
B( 8,3) = 0;
B( 9,3) = 0;
B(10,3) = 0;
B(11,3) = 0;
B(12,3) = 0;

% cwd (large wood)

B( 1,4) = 0;
B( 2,4) = 0;
B( 3,4) = 0;
B( 4,4) = 0;
B( 5,4) = 0;
B( 6,4) = 1;
B( 7,4) = 0;
B( 8,4) = 0;
B( 9,4) = 0;
B(10,4) = 0;
B(11,4) = 0;
B(12,4) = 0;

% cwd (coarse root)

B( 1,5) = 0;
B( 2,5) = 0;
B( 3,5) = 0;
B( 4,5) = 0;
B( 5,5) = 0;
B( 6,5) = 0;
B( 7,5) = 1;
B( 8,5) = 0;
B( 9,5) = 0;
B(10,5) = 0;
B(11,5) = 0;
B(12,5) = 0;
