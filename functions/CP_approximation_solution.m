function [Rho,info] = CP_approximation_solution(T,r,max_iter)
n = size(T);
d = length(n);
Rho = cell(d,1);

for mu = 1:d
    Rho{mu} = randn(n(mu),r);
end

L = cell(d,1);
L{1} = ones(1,r);
R = cell(d,1);
R{d} = ones(1,r);

for iter = 1:max_iter
    
    for mu = d-1:-1:1
        R{mu} = mode_1_kron(Rho{mu+1},R{mu+1});
    end
    
    for mu = 1:d
        LR = mode_1_kron(L{mu},R{mu});
        b = reshape(permute(T,[1:mu-1,mu+1:d,mu]),[prod([n(1:mu-1),n(mu+1:d)]),n(mu)]);
        Rho{mu} = (LR\b)';
        
        if mu < d
            L{mu+1} = mode_1_kron(L{mu},Rho{mu});
        end
    end
end
info = 0;
end