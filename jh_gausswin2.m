function [w] = jh_gausswin2(L,a)

% return a Gaussian window with alpha standard devisions.
%
% L .. number of points
% a .. number of std

N = L-1;
n = (0:N)'-N/2;
w = exp(-(1/2)*(a*n/(N/2)).^2);
end