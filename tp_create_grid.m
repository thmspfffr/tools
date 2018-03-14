  
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
  atlas.label_6mm = atlas.tissue_6mm(atlas.tissue_6mm>0)
  

  


