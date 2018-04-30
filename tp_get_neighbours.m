function neighbours = tp_get_neighbours(prj)
% put in the locations of sensors (either 2D or 3D) and the function
% returns the neighborhood structure based on delaunay triangulation
% the output format is NxN, where N is the number of channels/sensors/etc.

if ndims(prj) == 3
  
  tri = delaunay(prj(:,1), prj(:,2),prj(:,3));
  neighbours = zeros(size(prj,1),size(prj,1));
  
  % mark neighbours according to triangulation
  for i=1:size(tri, 1)
    
    neighbours(tri(i, 1), tri(i, 2)) = 1;
    neighbours(tri(i, 1), tri(i, 3)) = 1;
    neighbours(tri(i, 1), tri(i, 4)) = 1;
    
    neighbours(tri(i, 2), tri(i, 1)) = 1;
    neighbours(tri(i, 2), tri(i, 3)) = 1;
    neighbours(tri(i, 2), tri(i, 4)) = 1;
    
    neighbours(tri(i, 3), tri(i, 1)) = 1;
    neighbours(tri(i, 3), tri(i, 2)) = 1;
    neighbours(tri(i, 3), tri(i, 4)) = 1;
    
    neighbours(tri(i, 4), tri(i, 1)) = 1;
    neighbours(tri(i, 4), tri(i, 2)) = 1;
    neighbours(tri(i, 4), tri(i, 3)) = 1;
    
  end
  
else
  
  tri = delaunay(prj(:,1),prj(:,2));
  
  neighbours = zeros(size(prj,1),size(prj,1));
  % mark neighbours according to triangulation
  
  for i=1:size(tri, 1)
    
    neighbours(tri(i, 1), tri(i, 2)) = 1;
    neighbours(tri(i, 1), tri(i, 3)) = 1;
    
    neighbours(tri(i, 2), tri(i, 1)) = 1;
    neighbours(tri(i, 2), tri(i, 3)) = 1;
    
    neighbours(tri(i, 3), tri(i, 1)) = 1;
    neighbours(tri(i, 3), tri(i, 2)) = 1;
    
  end
end

% 

% save(sprintf('~/pconn/proc/preproc/neighb_sens274.mat'),'neighbours')

