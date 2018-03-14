
% ---------------------------
function X = jh_pinv(A,r)

if nargin<2
  r = 2;
end

[U,S,V] = svd(A,0);

s = diag(S);

s = diag(ones(r,1)./s(1:r));

X = V(:,1:r)*s*U(:,1:r)';


