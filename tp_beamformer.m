function [filt, pow] = tp_beamformer(Cf_real,lf,para)

% Compute a joint filter for all conditions

% input:
% Cf_real:  (real part of the complex) cross spectral density matrix
% lf:       leadfield
% para.reg: regularization parameter (typically para.reg = 0.05)

Cf_real = real(Cf_real);

% lambda = sum(diag(Cf_real))*para.reg; 
lambda = para.reg*sum(diag(Cf_real))/size(Cf_real,1); 

invCf_real = pinv(Cf_real + lambda * eye(size(Cf_real)));

n_voxel = size(lf,2);
% 
% full_filt = zeros(n_voxel,size(Cf_real,1),3);

for ilf=1:n_voxel 
  
%   disp(sprintf('Computing filter at location %d ...',ilf));
  
  filter = jh_pinv(squeeze(lf(:,ilf,:))' * invCf_real * squeeze(lf(:,ilf,:))) * squeeze(lf(:,ilf,:))' * invCf_real; % lf cell with leadfield, lf{i} : [chan x 3]

  cf = filter*Cf_real*filter';
  
  [u,s,~] = svd(cf);
   
  % extract power
  pow(ilf,1)=s(1,1);
  
  % projection to principal direction
  filt(ilf,:) = u(:,1)'*filter;
  
%   full_filt(ilf,:,:) = filter';

  
end

filt = filt';
