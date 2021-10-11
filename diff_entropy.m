function hc = diff_entropy(cov)
% compute differential entropy
% see e.g. mathworld.wolfram.com/DifferentialEntropy.html
% Takes covariance matrix as input, outputs diff entropy (hc)
% -------------

% cov(isnan(cov))=1;
s       = svd(cov);
logdetc = log(prod(s(s>0.00001)));
K       = sum(s>0.00001);
hc      = K/2*(1+log(2*pi))+1/2*logdetc;