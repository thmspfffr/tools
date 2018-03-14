function i = tp_matchaal(grid)

switch grid
  case 'coarse'
    load ~/pconn/matlab/aalmask_grid_coarse.mat
  case 'xcoarse'
    
  case 'cortex'
    
  case 'medium'
end

for imask = 1 : max(aalgrid.mask)
  
  i=find(aalgrid.mask==imask,1,'first');
  lab{imask}=cell2mat(aalgrid.labels(i));
  
end

f=fopen('~/Documents/MATLAB/aal_symm.nii.txt','rt');
a=textscan(f,'%d %s %d','headerlines',0);

for imask = 1 : 90
  
  if ~isempty(find(strcmp(a{2}(imask),lab)))
    i(imask)=find(strcmp(a{2}(imask),lab));
  else 
    i(imask)=nan;
  end
  
end

end