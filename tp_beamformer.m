function [filt, pow] = tp_beamformer(Cf_real,lf,para)

% Compute a joint filter for all conditions

% input:
% dat:    complex cross spectrum
% lf:     leadfield
% reg:    regularization parameter (default is set below)

% regularization parameter
lambda = sum(diag(Cf_real))*para.reg; 

invCf_real = pinv(Cf_real + lambda * eye(size(Cf_real)));

n_voxel = size(lf,2);

for ilf=1:n_voxel 
  
%   disp(sprintf('Computing filter at location %d ...',ilf));
  
  filter = jh_pinv(squeeze(lf(:,ilf,:))' * invCf_real * squeeze(lf(:,ilf,:))) * squeeze(lf(:,ilf,:))' * invCf_real; % lf cell with leadfield, lf{i} : [chan x 3]

  cf = filter*Cf_real*filter';
  
  [u,s,~] = svd(cf);
  
  % extract power
  pow(ilf)=s(1,1);
  % projection to principal direction
  filt(ilf,:) = u(:,1)'*filter; 
  
end

filt = filt';
