function [T,r] = TT_truncate_full_tensor(T,tol)

n = size(T);
d = length(n);
r = ones(1,d+1);

for mu = 1:d-1
    [U,S,V] = svd( reshape(T,[prod(n(1:mu)),prod(n(mu+1:d))]) ,0);
    r(mu+1) = sum(diag(S)>tol); % count number of entries larger than tol
    T = U(:,1:r(mu+1)) * S(1:r(mu+1),1:r(mu+1)) * V(:,1:r(mu+1))';
end

T = reshape(T,n);
end