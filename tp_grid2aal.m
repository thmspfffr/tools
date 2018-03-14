function aal = tp_grid2aal(val,para)
% converts input matrix to aal matrix
% value is connectivity matrix
% para needs
% - grid: (x)coarse/medium/cortex

switch para.grid
  
  case 'xcoarse'
    load ~/pconn/matlab/aalmask_grid_xcoarse
  case 'coarse'
    load ~/pconn/matlab/aalmask_grid_coarse
  case 'medium'
    load ~/pconn/matlab/aalmask_grid_medium
  case 'cortex'
    load('~/Documents/MATLAB/aalmask_grid_cortex3000.mat')
end

reg = unique(sort(aalgrid.mask(aalgrid.mask<=91)));
reg = reg(2:end);

if sum(size(squeeze(val))>1)==2 && size(val,1)==size(val,2)
  
  for ireg = 1 : length(reg)

    idx1 = aalgrid.mask==reg(ireg);

    for jreg = 1 : length(reg)

      idx2 = aalgrid.mask==reg(jreg);

      if reg(ireg) == reg(jreg)
        aal(ireg,jreg) = NaN;
      else
        aal(ireg,jreg) = nanmean(nanmean(val(idx1,idx2),1),2);
      end

    end
  end
elseif sum(size(squeeze(val))>1) == 3 && size(val,1)==size(val,2)
  
  for ireg = 1 : length(reg)

    idx1 = aalgrid.mask==reg(ireg);

    for jreg = 1 : length(reg)

      idx2 = aalgrid.mask==reg(jreg);

      if reg(ireg) == reg(jreg)
        aal(ireg,jreg,:) = nan(2,1);
      else
        aal(ireg,jreg,:) = nanmean(nanmean(val(idx1,idx2,:),1),2);
      end

    end
  end
elseif sum(size(squeeze(val))>1) == 2 && size(val,1)~=size(val,2)
  
  for ireg = 1 : length(reg)
    
    idx1 = aalgrid.mask==reg(ireg);
    
    aal(:,ireg) = nanmean(val(:,idx1),2);

  end
  
else
  error('Data are neither 2D or 3D! Check data format!')
end

% aal.regnr = reg;
% aal.lab = aalgrid.labels(reg);
