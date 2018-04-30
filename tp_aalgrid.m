function [aal] = tp_aalgrid
%function mask = grid2aalmask(atlas, grid)

% creates grid based on AAL regions, outputs AAL positions for creating
% leadfield in pconn_src_sa.m

fprintf('Calculating atlas region mask. This might take a while....\n');

aal = ft_read_atlas('~/Documents/MATLAB/fieldtrip-20160919/template/atlas/aal/ROI_MNI_V4.nii');
% get coordinates from voxel indices
ac=zeros(aal.dim(1)*aal.dim(2)*aal.dim(3),3);

ii=0;
for ix=1:aal.dim(1)
  for iy=1:aal.dim(2)
    for iz=1:aal.dim(3)
      ii=ii+1;
      r=[ix,iy,iz,1]';
      rc=aal.transform*r;
      ac(ii,:)=rc(1:3)';
      aal.tissue_2mm(ii) = aal.tissue(ix,iy,iz);
    end
  end
end

aal.grid_2mm = ac(aal.tissue_2mm>0&aal.tissue_2mm<=91,:);
% 

% subsample
ii=0;
for ix=1:3:aal.dim(1)
  for iy=1:3:aal.dim(2)
    for iz=1:3:aal.dim(3)
      ii=ii+1;
      r=[ix,iy,iz,1]';
      rc=aal.transform*r;
      bc(ii,:)=rc(1:3)';
      aal.tissue_6mm(ii) = aal.tissue(ix,iy,iz);
    end
  end
end

aal.grid_6mm = bc(aal.tissue_6mm>0&aal.tissue_6mm<=91,:);

clear bc

ii=0;
for ix=1:2:aal.dim(1)
  for iy=1:2:aal.dim(2)
    for iz=1:2:aal.dim(3)
      ii=ii+1;
      r=[ix,iy,iz,1]';
      rc=aal.transform*r;
      bc(ii,:)=rc(1:3)';
      aal.tissue_4mm(ii) = aal.tissue(ix,iy,iz);
    end
  end
end

aal.grid_4mm = bc(aal.tissue_4mm>0&aal.tissue_4mm<=91,:);

% compute center of mass of each AAL region
% ----------------
for i = 1 : 91
  aal.centroids(i,:) = mean(aal.grid_2mm(aal.tissue_2mm(aal.tissue_2mm>0&aal.tissue_2mm<=91)==i,:));
end
