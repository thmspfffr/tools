  
function [m758] = tp_create_grid(para)

fprintf('Calculating atlas region mask...\n');

if strcmp(para,'m758')
  load ~/Documents/MATLAB/m758_atlas.mat
  atlas = m758;
elseif strcmp(para,'vtpm')
  load ~/Documents/MATLAB/m758_atlas.mat
  atlas = vtpm;
end
  % get coordinates from voxel indices
  ac=zeros(atlas.dim(1)*atlas.dim(2)*atlas.dim(3),3);

  clear bc

  ii=0;
  for ix=1:2:atlas.dim(1)
      for iy=1:2:atlas.dim(2)
          for iz=1:2:atlas.dim(3)
              ii=ii+1;
              r=[ix,iy,iz,1]';
              rc=atlas.transform*r;
              bc(ii,:)=rc(1:3)'; 
              atlas.tissue_4mm(ii) = atlas.anatomy(ix,iy,iz);
          end
      end
  end

  atlas.grid_4mm = bc(atlas.tissue_4mm>0,:);

  clear bc

  ii=0;
  for ix=1:3:atlas.dim(1)
      for iy=1:3:atlas.dim(2)
          for iz=1:3:atlas.dim(3)
              ii=ii+1;
              r=[ix,iy,iz,1]';
              rc=atlas.transform*r;
              bc(ii,:)=rc(1:3)'; 
              atlas.tissue_6mm(ii) = atlas.anatomy(ix,iy,iz);
          end
      end
  end

  atlas.grid_6mm = bc(atlas.tissue_6mm>0,:);
  

  


