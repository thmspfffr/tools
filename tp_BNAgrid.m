function [BNA] = tp_BNAgrid
%function mask = grid2aalmask(atlas, grid)

% creates grid based on AAL regions, outputs AAL positions for creating
% leadfield in pconn_src_sa.m

restoredefaultpath
fprintf('Calculating atlas region mask. This might take a while....\n');
addpath ~/Documents/MATLAB/fieldtrip-20190224/
ft_defaults
BNA = ft_read_atlas('~/Documents/MATLAB/fieldtrip-20190224/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii');

% get coordinates from voxel indices
ac=zeros(BNA.dim(1)*BNA.dim(2)*BNA.dim(3),3);

% subsample: 4 = 5mm (5*1.25mm)
subsamp = 4;

ii=0;
for ix=1:subsamp:BNA.dim(1)
  for iy=1:subsamp:BNA.dim(2)
    for iz=1:subsamp:BNA.dim(3)
      ii=ii+1;
      r=[ix,iy,iz,1]';
      rc=BNA.transform*r;
      bc(ii,:)=rc(1:3)';
      BNA.tissue_5mm(ii) = BNA.tissue(ix,iy,iz);
    end
  end
end

BNA.grid_5mm = bc(BNA.tissue_5mm>0&BNA.tissue_5mm<=length(BNA.tissuelabel),:);
BNA.tissue_5mm =  BNA.tissue_5mm(BNA.tissue_5mm>0&BNA.tissue_5mm<=length(BNA.tissuelabel))
% compute center of mass of each AAL region
% ----------------
for i = 1 : length(BNA.tissuelabel)
  BNA.centroids(i,:) = mean(BNA.grid_5mm(BNA.tissue_5mm(BNA.tissue_5mm>0&BNA.tissue_5mm<=length(BNA.tissuelabel))==i,:));
end

