function geneidx = tp_aalgenemask(nreg)

addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')
ft_defaults
aal = ft_read_atlas('~/Documents/MATLAB/fieldtrip-20160919/template/atlas/aal/ROI_MNI_V4.nii');
load ~/pmod/matlab/gene_values.mat

geneidx = [];

for i = 1:nreg
  idx = find(aal.tissue==i);
  
  for ii = 1 : length(idx)
    [x,y,z]=ind2sub(aal.dim,idx(ii));
    
    tmp = find(sum(locs==[x y z],2)==3);
    if ~isempty(tmp)
      geneidx = [geneidx(:); tmp];
    end
    
  end
end
