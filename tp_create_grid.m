
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
  %   stepsize = [4 6];
  stepsize = [2 3];
else
  stepsize = [2 3];
end
% get coordinates from voxel indices
ac=zeros(atlas.dim(1)*atlas.dim(2)*atlas.dim(3),3);

clear bc

% label correctly
for i = 1 : 50
  if i < 26; prefix = 'l_'; else prefix = 'r_'; end
  atlas.tissuelabel{i} = [prefix atlas.tissuelabel{i}];
end

% ----------------------------------------
% stepsize(1): 4 mm?
% ----------------------------------------

ii=0; bc = zeros(prod(floor(atlas.dim./stepsize(1))),3);
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

% get left/right
% ----------------------------------------
for i = 1 : size(atlas.grid_4mm,1)
  if atlas.grid_4mm(i,1)<0
    add_const = 0;
  else
    add_const = 25;
  end
  atlas.new_label_4mm(i) = atlas.label_4mm(i)+add_const;
end
atlas.label_4mm = atlas.new_label_4mm;
atlas = rmfield(atlas,{'new_label_4mm'});

% delete empty stuff
% ----------------------------------------
l_idx = unique(atlas.label_4mm(atlas.label_4mm<26));
r_idx = unique(atlas.label_4mm(atlas.label_4mm>25))-25;
lr = intersect(l_idx,r_idx);
idx = [lr lr+25];
% atlas.tissuelabel_4mm = [atlas.tissuelabel(idx(1:length(idx)/2)) atlas.tissuelabel(idx(length(idx)/2+1:end))];
atlas.tissuelabel_4mm = atlas.tissuelabel(idx);

for i = 1 : 50
  old_idx{i} = find(atlas.label_4mm==i);
end

lab = zeros(1,size(atlas.label_4mm,2));

for i = 1 : length(idx)
  lab(old_idx{idx(i)})=i;
end

atlas.label_4mm = lab; 
  
% figure; set(gcf,'color','w')

clear bc

% ----------------------------------------
% stepsize(2): 6 mm?
% ----------------------------------------

ii=0; bc = zeros(prod(floor(atlas.dim./stepsize(2))),3);
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

% get left/right
% ----------------------------------------
for i = 1 : size(atlas.grid_6mm,1)
  if atlas.grid_6mm(i,1)<0
    add_const = 0;
  else
    add_const = 25;
  end
  atlas.new_label_6mm(i) = atlas.label_6mm(i)+add_const;
end
atlas.label_6mm = atlas.new_label_6mm;
atlas = rmfield(atlas,{'new_label_6mm'});

% delete empty stuff
% ----------------------------------------
l_idx = unique(atlas.label_6mm(atlas.label_6mm<26));
r_idx = unique(atlas.label_6mm(atlas.label_6mm>25))-25;
lr = intersect(l_idx,r_idx);
not = setdiff([l_idx r_idx],lr);
idx = [lr lr+25];

atlas.tissuelabel_6mm = atlas.tissuelabel(idx);

for i = 1 : 50
  old_idx{i} = find(atlas.label_6mm==i);
end

lab = zeros(1,size(atlas.label_6mm,2));

for i = 1 : length(idx)
  lab(old_idx{idx(i)})=i;
end

atlas.label_6mm = lab; 
clear bc
