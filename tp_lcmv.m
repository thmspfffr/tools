function [lcmv, A]=tp_lcmv(L,C,alpha0)
% Calculates spacial filters, and  power over all grid points using lcmv
% Input:
% L: Lead field matrix, NxMx3 matrix for N channels, M voxels and 3 dipole
%    directions
% C:  NxN matrix for N channels covariance matrix or real part of cross-spectrum
%    (The program uses only the real of C)
% alpha0: (relative) regularization parameter. In the algorthm C+alpha*eye(N) is
%        inverted with alpha=alpha0*trace(C)/N
%
% Output
% A  :3-dimensional filter
% lcmv.filt : 1-dimensional filter along direction with strongest power
% lcmv.pow: Mx1 vector for M voxels, po(i) is the power at the i.th voxel along
%          strongest direction
% lcmv.noise: noise projections for computing NAI: (par.pow./par.noise)

C=real(C);

if nargin<3
    alpha0=.05;
end

alpha=alpha0*trace(C)/length(C);
 
[nchan nchan]=size(C);
[nchan ns ndum]=size(L);
Cr=C+alpha*eye(nchan);


noise = svd(C);
noise = noise(end);
% estimated noise floor is equal to or higher than lambda
noise = max(noise, alpha);
    
Crinv=inv(Cr);

A=zeros(nchan,ns,3);
for i=1:ns
    Lloc=squeeze(L(:,i,:));
    A(:,i,:)=reshape((inv((Lloc'*Crinv*Lloc))*Lloc'*Crinv)',nchan,3);
end


po=zeros(ns,1);
lcmv.filt=zeros(nchan,ns);
for i=1:ns
  i
    Aloc                = transpose(squeeze(A(:,i,:)));
    Ploc                = Aloc * C * Aloc';
    [u s v]             = svd(Ploc);
    lcmv.pow(i,1)        = s(1,1);
    lcmv.filt(:,i)     	= Aloc'*u(:,1);
    s                   = svd(lcmv.filt(:,i) * lcmv.filt(:,i)'); 
    x                   = s(1);
    lcmv.noise(i,1)      = noise * x;
    lcmv.noisetrace(i,1) = noise*trace(lcmv.filt(:,i)*lcmv.filt(:,i)');
end




return;
