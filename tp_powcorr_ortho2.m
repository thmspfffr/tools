function [c] = tp_powcorr_ortho(data,para,sa)

% computes hipps stuff

if ~exist('para','var')
  para.segave = 1;
end

epleng    = para.epleng;
epshift   = para.epshift;
segleng   = para.segleng;
fsample   = para.fsample;
segshift  = para.segshift;
f         = para.foi;

scale = sqrt(mean(mean(data.^2)));
data  = data/scale;

% [~, ns1] = size(A1);
% [~, ns2] = size(A2);

[n nchan] = size(data);

switch para.wavelet
  case 'hanning'
    nn    = (1:segleng)'-segleng/2;
    mywin = hanning(segleng);
    s1 = cos(nn*f*2*pi/fsample).*mywin;
    s2 = sin(nn*f*2*pi/fsample).*mywin;
    ss = s1-sqrt(-1)*s2;
    ss = repmat(ss,1,nchan);
  case 'ft'
    para.fsample = fsample;
    para.freqoi = f;
    para.gwidth = 3;
    para.width  = 7;
    
    w = tp_mkwavelet(para);
    ss = repmat(w,1,nchan);
    segleng = length(w);
    segshift = floor(segleng/2);
    
end

% 
% if nargin<8;
%   A2 = A1;
% end

nep  = floor((n-epleng)/epshift+1);
nseg = floor((epleng-segleng)/segshift+1);

for iep = 1 : nep
  
  dloc1 = data((iep-1)*epshift+1:(iep-1)*epshift+epleng,:);
  
  for iseg = 1 : nseg
    
    fprintf('Seg %d/%d ...\n',iseg,nseg)
    
%     dloc2 = dloc1((iseg-1)*segshift+1:(iseg-1)*segshift+segleng,:);
    dloc2 = fft(dloc1((iseg-1)*segshift+1:(iseg-1)*segshift+segleng,:));
    wfft  = fft(w);
    
    dum = fftshift(ifft(dloc2 .* repmat(wfft,[1 nchan]), [], 2),2);
    dum = dum .* sqrt(2 ./ fsample);
    dum = mean(dum);

%     dataf = sum(dloc2.*ss);
    
%     mom(:,iseg) = dataf*A1;
    mom(:,iseg) = dum;
  end
  
  filt = pconn_beamformer(mom,sa.sa.L_aal);
  
  mom = filt'*mom;
  
  % this one computes orthopowcorrs based on fieldtrip code
  c(:,:,iep) =  compute_orthopowcorr(mom); clear mom
  
end
% compute power corr


function c = compute_orthopowcorr(mom,varargin)

refindx = 'all'; clear tapvec
tapvec  = ones(1,size(mom,2));

if strcmp(refindx, 'all')
  refindx = 1:size(mom,1);
end

cmomnorm = conj(mom./abs(mom)); % only need to do conj() once

n        = size(mom,1);
ntap     = 1;
if ~all(tapvec==ntap)
  error('unequal number of tapers per observation is not yet supported');
end

ix = zeros(sum(tapvec),1);
jx = ix;
sx = ix;
for k = 1:numel(tapvec)
  indx = (k-1)*ntap+(1:ntap);
  ix(indx) = indx;
  jx(indx) = k;
  sx(indx) = 1./ntap;
end
tra = sparse(ix,jx,sx,sum(tapvec),numel(tapvec));

powmom = (abs(mom).^2)*tra; % need only once
powmom = standardise(log10(powmom), 2);

c = zeros(n, numel(refindx));%;*2);
N = ones(n,1);
%warning off;
for k = 1:numel(refindx)
  indx     = refindx(k);
  ref      = mom(indx,:);
  crefnorm = conj(ref./abs(ref));
  
  pow2 = (abs(imag(ref(N,:).*cmomnorm)).^2)*tra;
  pow2 = standardise(log10(pow2), 2);
  c1   = mean(powmom.*pow2, 2);
  pow1 = (abs(imag(mom.*crefnorm(N,:))).^2)*tra;
  pow1 = standardise(log10(pow1), 2);
  
  pow2 = (abs(ref).^2)*tra;
  pow2 = standardise(log10(pow2), 2);
  pow2 = repmat(pow2, [n 1]);
  c2   = mean(pow1.*pow2, 2);
  
  c(:,k) = (c1+c2)./2;
  %c(:,k+numel(refindx)) = c2;
end


