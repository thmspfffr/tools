function [filt, pow,cf] = tp_beamformer(Cf_real,lf,para)
% tp_beamformer yields a spatial filter according to LCMV (van Veen et al.,
% 1997).
% ----------
% INPUT
% ----------
% Cf_real:  (real part of the complex) cross spectral density matrix
% lf:       leadfield
% para.reg: regularization parameter (typically para.reg = 0.05)
% ----------
% OUTPUT
% ----------
% filt: 1D spatial filter
% pow: source level power estimates

Cf_real     = real(Cf_real);
lambda      = para.reg*sum(diag(Cf_real))/size(Cf_real,1); 
invCf_real  = pinv(Cf_real + lambda * eye(size(Cf_real)));

n_voxel = size(lf,2);

filt = zeros(n_voxel,size(Cf_real,1));
for ilf=1:n_voxel 
    
  filter  = jh_pinv(squeeze(lf(:,ilf,:))' * invCf_real * squeeze(lf(:,ilf,:))) * squeeze(lf(:,ilf,:))' * invCf_real; % lf cell with leadfield, lf{i} : [chan x 3]
  cf      = filter*Cf_real*filter';
  [u,s,~] = svd(cf);
  % projection to principal direction
  filt(ilf,:) = u(:,1)'*filter;
  % extract source level power estimates
  pow(ilf,1)  = s(1,1);
  
end

filt = filt';
