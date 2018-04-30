%% AAL GET REGIONS NAMES

function [transdat, trans, names] = tp_match_aal(para,varargin)

clear names trans reg

switch para.grid
  case 'cortex'
    load aalmask_grid_cortex3000
  case 'coarse'
    load aalmask_grid_coarse
  case 'medium'
    load aalmask_grid_medium
    case 'xcoarse'
    load aalmask_grid_xcoarse
end

if ~isfield(para,'N')
  pars.N = 91;
end

reg = unique(sort(aalgrid.mask(aalgrid.mask<=para.N)));

if reg(1) == 0
  reg = reg(2:end);
end

for ireg = 1 : length(reg)
  
  iireg = reg(ireg);
  
  names(iireg) = unique(aalgrid.labels(aalgrid.mask==iireg));
  
end

% MATCH WITH ATLAS FROM BCN

f=fopen('~/Documents/MATLAB/aal_symm.nii.txt','rt');
a=textscan(f,'%d %s %d','headerlines',0);

for imask = 1 : length(reg)
  
  iireg = reg(imask);
  
  if ~isempty(find(strcmp(a{2},names(iireg))))
    trans(imask,:)=[imask find(strcmp(a{2},names(iireg)))];
  else
    trans(imask,:)=[imask nan];
  end
  
end

if exist('varargin','var')
  
  transdat = nan(90,90,size(varargin{1},3));
  
  dat                        = varargin{1};
  dat(isnan(trans(:,2)),:,:) = [];
  dat(:,isnan(trans(:,2)),:) = [];
  trans(isnan(trans(:,2)),:) = [];
  
  for i = 1 : size(trans)
    for j = 1 : size(trans)
      
      transdat(trans(i,2),trans(j,2),:) = dat(trans(i,1),trans(j,1),:);
    
    end
  end
end




