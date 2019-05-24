function [r,p] = tp_corr(x,y,dim)
% COMPUTE PEARSON CORRELATION FOR MULTI-DIM ARRAYS
% ---------------------------------
% The function computes pearson correlation coefficient and p-values
% (derived from a t-distribution) for multi-dimensional arrays. 
% Inputs:
% - x,y: data to correlate
% - dim: dimension across which data should be correlated
% -----------
% Example:
% Let x and y be a 40x40x30 array, where the correlation should be computed
% across the third dimension. 
% The function should be called: tp_corr(x,y,3);
% -----------

if dim ~= 3
  error('Not tested for non 3D arrays')
end
  
az = bsxfun(@minus, x, nanmean(x,dim));
bz = bsxfun(@minus, y, nanmean(y,dim));
  
a2 = az .^ 2;
b2 = bz .^ 2;
ab = az .* bz;

r = sum(ab,3) ./ sqrt(sum(a2,dim) .* sum(b2,dim));
t = r.* sqrt ( (28-2 ) ./ 1 - r.^2);

p=2*tcdf(-abs(t),27);
  