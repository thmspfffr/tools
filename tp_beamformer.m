function [filt pow NAI] = pconn_beamformer(dat,lf,para)

% Compute a joint filter for all conditions

% input:
% dat:    complex cross spectrum
% lf:     leadfield
% reg:    regularization parameter (default is set below)

if para.iscs
  Cf_real = real(dat);
else
 % dat is complex freq-specific data [chan x samples]
  Cf_real = real(dat*dat'/size(dat,2));
end
% regularization parameter
lambda = sum(diag(Cf_real))*para.reg; 

invCf_real = pinv(Cf_real + lambda * eye(size(Cf_real)));

clear filt

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
%   NAI(ilf) = trace((squeeze(lf(:,ilf,:))'*invCf_real* squeeze(lf(:,ilf,:))).^(-1)) / trace(squeeze(lf(:,ilf,:))'*para.noisecov*squeeze(lf(:,ilf,:)));

end

filt = filt';
