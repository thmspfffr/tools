  
function atlas = tp_create_grid(para)

fprintf('Calculating atlas region mask...\n');

if strcmp(para,'m758')
  load ~/Documents/MATLAB/m758_atlas.mat
  atlas = m758;
elseif strcmp(para,'vtpm')
  load ~/Documents/MATLAB/fieldtrip-20160919/template/atlas/vtpm/vtpm.mat
  atlas = vtpm;
  atlas.anatomy = atlas.tissue;
end

if prod(size(atlas.anatomy)) == 7221032
  stepsize = [4 6];
else
  stepsize = [2 3];
end
  % get coordinates from voxel indices
  ac=zeros(atlas.dim(1)*atlas.dim(2)*atlas.dim(3),3);

  clear bc
  
  ii=0;
  for ix=1:stepsize(1):atlas.dim(1)
      for iy=1:stepsize(1):atlas.dim(2)
          for iz=1:stepsize(1):atlas.dim(3)
              ii=ii+1;
              r=[ix,iy,iz,1]';
              rc=atlas.transform*r;
              bc(ii,:)=rc(1:3)'; 
              atlas.tissue_4mm(ii) = atlas.anatomy(ix,iy,iz);
          end
      end
  end

  atlas.grid_4mm = bc(atlas.tissue_4mm>0,:);
  atlas.label_4mm = atlas.tissue_4mm(atlas.tissue_4mm>0);
  clear bc

  ii=0;
  for ix=1:stepsize(2):atlas.dim(1)
      for iy=1:stepsize(2):atlas.dim(2)
          for iz=1:stepsize(2):atlas.dim(3)
              ii=ii+1;
              r=[ix,iy,iz,1]';
              rc=atlas.transform*r;
              bc(ii,:)=rc(1:3)'; 
              atlas.tissue_6mm(ii) = atlas.anatomy(ix,iy,iz);
          end
      end
  end

  atlas.grid_6mm = bc(atlas.tissue_6mm>0,:);
  atlas.label_6mm = atlas.tissue_6mm(atlas.tissue_6mm>0);
  
%% IMPLEMENT AAL GRID HERE
% -----------------------------------
  
% 
% aal = ft_read_atlas('~/Documents/MATLAB/fieldtrip-20160919/template/atlas/aal/ROI_MNI_V4.nii');
% % get coordinates from voxel indices
% ac=zeros(aal.dim(1)*aal.dim(2)*aal.dim(3),3);
% 
% ii=0;
% for ix=1:aal.dim(1)
%   for iy=1:aal.dim(2)
%     for iz=1:aal.dim(3)
%       ii=ii+1;
%       r=[ix,iy,iz,1]';
%       rc=aal.transform*r;
%       ac(ii,:)=rc(1:3)';
%       aal.tissue_2mm(ii) = aal.tissue(ix,iy,iz);
%     end
%   end
% end
% 
% % subsample
% ii=0;
% for ix=1:3:aal.dim(1)
%   for iy=1:3:aal.dim(2)
%     for iz=1:3:aal.dim(3)
%       ii=ii+1;
%       r=[ix,iy,iz,1]';
%       rc=aal.transform*r;
%       bc(ii,:)=rc(1:3)';
%       aal.tissue_6mm(ii) = aal.tissue(ix,iy,iz);
%     end
%   end
% end
% 
% aal.grid_6mm = bc(aal.tissue_6mm>0&aal.tissue_6mm<=91,:);
% 
% clear bc
% 
% ii=0;
% for ix=1:2:aal.dim(1)
%   for iy=1:2:aal.dim(2)
%     for iz=1:2:aal.dim(3)
%       ii=ii+1;
%       r=[ix,iy,iz,1]';
%       rc=aal.transform*r;
%       bc(ii,:)=rc(1:3)';
%       aal.tissue_4mm(ii) = aal.tissue(ix,iy,iz);
%     end
%   end
% end
% 
% aal.grid_4mm = bc(aal.tissue_4mm>0&aal.tissue_4mm<=91,:);
% 
% 
% % compute center of mass of each AAL region
% 
% for i = 1 : size(aal.tissuelabel,2)
%   
%   aal.centroids(i,:) = mean(ac(aal.tissue==i,:));
%   
% end


