function U = decorrgs(V,varargin)    
% U = decorrgs(V,dim)  
% Orthonormalise a set of vectors via the Gram-Schmidt process 


if nargin == 1
    dim = 1;
else
    dim = varargin{1};
end


if dim == 1
    V = V';
end


n = size(V,1);
k = size(V,2);
U = zeros(n,k);
U(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
for i = 2:k
    U(:,i) = V(:,i);
    for j = 1:i-1
        U(:,i) = U(:,i) - ( U(:,i)'*U(:,j) )/( U(:,j)'*U(:,j) )*U(:,j);
    end
    U(:,i) = U(:,i)/sqrt(U(:,i)'*U(:,i));
end

if dim == 1
    U = U';
end