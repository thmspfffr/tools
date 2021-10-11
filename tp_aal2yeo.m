function [lab,yeo1mm] = tp_grid2yeo(grid)

yeo1mm = ft_read_atlas('~/Documents/MATLAB/Yeo2011_7Networks_MNI152_FreeSurferConformed1mm.nii.gz')

ii=0;

yeo1mm.brick0_2mm = zeros(128*128*128,1);
bc = zeros(128*128*128,3);

for ix=2:2:yeo1mm.dim(1)
  for iy=2:2:yeo1mm.dim(2)
    for iz=2:2:yeo1mm.dim(3)
      
      ii=ii+1;
      r=[ix,iy,iz,1]';
      rc=yeo1mm.transform*r;
      bc(ii,:)=rc(1:3)';
      yeo1mm.brick0_2mm(ii) = yeo1mm.brick0(ix,iy,iz);
      
    end
  end 
end

yeo1mm.grid_2mm = bc(yeo1mm.brick0_2mm>0,:);
yeo1mm.label_2mm = yeo1mm.brick0_2mm(yeo1mm.brick0_2mm>0);
% 
% grid=ft_read_atlas('~/Documents/MATLAB/fieldtrip-20160919/template/atlas/grid/ROI_MNI_V4.nii')
% 
% ii = 0;
% grid.brick0_2mm = zeros(grid.dim(1)*grid.dim(2)*grid.dim(3),1);
% bc = zeros(grid.dim(1)*grid.dim(2)*grid.dim(3),3);
% 
% for ix=1:grid.dim(1)
%   for iy=1:grid.dim(2)
%     for iz=1:grid.dim(3)
% 
%       ii=ii+1;
%       r=[ix,iy,iz,1]';
%       rc=grid.transform*r;
%       bc(ii,:)=rc(1:3)';
%       grid.brick0_2mm(ii) = grid.tissue(ix,iy,iz);
%       
%     end
%   end
% end
% 
% grid.grid_2mm = bc(grid.brick0_2mm>0,:);
% grid.label_2mm = grid.brick0_2mm(grid.brick0_2mm>0);


% for i = 1 : 90
%   i
%   idx = find(grid.label_2mm == i);
%   iii = 0;
%   for ii = 1 : length(idx)
%     if ~isempty(find(sum((grid.grid_2mm(idx(ii),:)-yeo1mm.grid_2mm)==0,2)==3))
%       iii = iii + 1;
%       idx_grid{i}(iii)=find(sum((grid.grid_2mm(idx(ii),:)-yeo1mm.grid_2mm)==0,2)==3);
%       sum_grid(i,iii)=sum(yeo1mm.label_2mm(idx_grid{i})==ii);
%     end
%   end
% end

for i = 1 : size(grid,1)
    idx_grid = 10*grid(i,:);
%   for ii = 1 : 7
%     idx_yeo = find(yeo1mm.label_2mm == ii);
    dist = sqrt(sum((idx_grid-yeo1mm.grid_2mm).^2,2));
    [~,min_dist_idx]=min(dist);
    lab(i) =  yeo1mm.label_2mm(min_dist_idx);
end
    
  
  

% 
% 
% for i = 1 : size(grid,1)
%     idx_grid = grid(i,:);
% %   for ii = 1 : 7
% %     idx_yeo = find(yeo1mm.label_2mm == ii);
%     j = 0;
%     for ii = 1 : length(idx_grid)
%       tmp_idx = find(sum((grid.grid_2mm(idx_grid(ii,:),:)-yeo1mm.grid_2mm)==0,2)==3);
%       if ~isempty(tmp_idx)
%         j = j + 1;
%         trans(j) = yeo1mm.label_2mm(tmp_idx);
%       end
%     end
%     
%     for irsn = 1 : 7
%     grid2yeo(i,irsn) = sum(trans==irsn);
%     end
% end
%     
% grid2yeo=grid2yeo./max(grid2yeo)
%   
% 
%   