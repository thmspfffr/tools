%% AAL GET REGIONS NAMES
% --------------------------
% Function transforms AAL order from "Hamburg layout" (which corresponds to
% the default layout (LR LR LR... order) to the "Barcelona layout" (which
% corresponds to (LLL... RRR order). If para.N is smaller than 90, the
% script will additionally eliminate some "subcortical" (and other
% unwanted) regions that are defined in exclude_bcn and exclude_hh.
% Definition of the grid doesn't really matter, its only purpuse is to
% extract the order of the HH atlas.
% -------------------------
% [transdat, trans, names] = tp_match_aal(para,varargin)
% para.grid = medium (default)
% para.N = 91 (default)
% para.transfer = 'to_hh' or 'to_bcn'
% para.dim = 2 (default), also words for vectors (e.g. power estimates)

function [transdat, trans, translab] = tp_match_aal(para,varargin)

% addpath ~/Documents/MATLAB/fieldtrip-20160919/; ft_defaults
aal = ft_read_atlas('~/Documents/MATLAB/fieldtrip-20160919/template/atlas/aal/ROI_MNI_V4.nii');

% exclude those regions in BCN definition
exclude_bcn = [11 15 16 17 18 19 20 21 36 37 38 39 52 53 54 55 70 71 72 73 74 75 76 80];
% exclude those regions in HH definition
exclude_hh  = [21 22 29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];
% ------

clear names trans reg
addpath ~/pconn/matlab

% if ~isfield(para,'grid')
%   para.grid = 'medium';
%   load aalmask_grid_medium
% else
%   switch para.grid
%     case 'cortex'
%       load aalmask_grid_cortex3000
%     case 'coarse'
%       load aalmask_grid_coarse
%     case 'medium'
%       load aalmask_grid_medium
%     case 'xcoarse'
%       load aalmask_grid_xcoarse
%   end
% end

if ~isfield(para,'N')
  para.N = 91;
end
if ~isfield(para,'dim')
  para.dim = 2;
end

if strcmp(para.transfer, 'to_hh')
  
  % MATCH WITH ATLAS FROM BCN
  
  f=fopen('~/Documents/MATLAB/aal_symm.nii.txt','rt');
  a=textscan(f,'%d %s %d','headerlines',0);
  
  for imask = 1 : para.N
    
%     iireg = reg(imask);
    
    if ~isempty(find(strcmp(a{2},aal.tissuelabel(imask))))
      trans(imask,:)=[imask find(strcmp(a{2},aal.tissuelabel(imask)))];
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
      translab{i} = a{2}{trans(i,2)};
      for j = 1 : size(trans)
        
%         transdat(trans(i,2),trans(j,2),:) = dat(trans(i,1),trans(j,1),:);
          transdat(i,j,:) = dat(trans(i,2),trans(j,2),:);
         
 
      end
    end
  end
  
  if para.dim == 1
    
    transdat = transdat(:,1);
    
  end
  
elseif strcmp(para.transfer,'to_bcn')
  
  f=fopen('~/Documents/MATLAB/aal_symm.nii.txt','rt');
  a=textscan(f,'%d %s %d','headerlines',0);
  names=a{2}';
  
  % MATCH WITH ATLAS FROM BCN
  for imask = 1 : length(names)
        
    if ~isempty(find(strcmp(aal.tissuelabel,names(imask))))
      trans(imask,:)=[imask find(strcmp(aal.tissuelabel,names(imask)))];
    else
      trans(imask,:)=[imask nan];
    end
    
  end
  transdat = nan(90,90,size(varargin{1},3));
    
    dat                        = varargin{1};
    dat(isnan(trans(:,2)),:,:) = [];
    dat(:,isnan(trans(:,2)),:) = [];
    trans(isnan(trans(:,2)),:) = [];
    
    for i = 1 : size(trans,1)
      translab{i} = aal.tissuelabel{trans(i,2)};
      for j = 1 : size(trans,1)
        transdat(i,j,:) = dat(trans(i,2),trans(j,2),:);
        
      end
    end
  
  if para.dim == 1
    
    transdat = transdat(:,1);
    
  end
end

idx = 1 : 90;

if para.N < 90 && para.dim == 2
  warning('Excluding regions...')
  if strcmp(para.transfer,'to_bcn')
    % this means, atlas is now in BCN format
    transdat = transdat(~ismember(idx,exclude_bcn),~ismember(idx,exclude_bcn));
  elseif strcmp(para.transfer,'to_hh')
    % this means, atlas is now in HH format (default format)
    transdat = transdat(~ismember(idx,exclude_hh),~ismember(idx,exclude_hh));
  end
elseif para.N < 90 && para.dim == 1
  warning('Excluding regions...')
  if strcmp(para.transfer,'to_bcn')
    % this means, atlas is now in BCN format
    transdat = transdat(~ismember(idx,exclude_bcn));
  elseif strcmp(para.transfer,'to_hh')
    % this means, atlas is now in HH format (default format)
    transdat = transdat(~ismember(idx,exclude_hh));
  end
end
  


