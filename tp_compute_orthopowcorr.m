function c = tp_compute_orthopowcorr(mom,varargin)

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
  tic
%   k
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
  toc
end