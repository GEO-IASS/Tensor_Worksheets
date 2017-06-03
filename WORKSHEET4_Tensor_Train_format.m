%% Tensor Train format 
% _Written by Sebastian Kraemer, IGPM at RWTH Aachen University_
%
% _The TT decomposition was first introduced to mathematics by I.V.
% Oseledets (2011)_
%
% The Tensor Train (short: TT) format may appear complicated, but in many
% aspects it is one of the most convenient formats.
%
% $$ $$ Similar to the Tucker format, the TT-ranks $r \in R^{d-1}$ are
% easily calculated as $$ $$
%
% $$ r_{\mu}(T) = \mathrm{TT-rank}_{\mu}(T) = \mathrm{rank}(T^{(1,\ldots,\mu)}) $$
%
% $$ $$ where the matrix $$ T^{(1,\ldots,\mu)} $$ is defined via (again using Matlab code)
% $$ $$
%
% $$ T^{(1,\ldots,\mu)} :=
% {\tt
% reshape(T,[1:mu,mu+1:d])}
% $$
%
% $$ $$ For simplicity, one sets $r_0 = r_d = 1$. It then follows that there exist cores $$ $$
%
% $$ G_\mu: \{1,\ldots,n_\mu\} \mapsto R^{r_{\mu-1} \times r_{\mu}},\ \mu = 1,\ldots,d $$ such that $$ T_I = G_1(i_1) \cdot \ldots \cdot G_d(i_d), \quad \forall I = (i_1,\ldots,i_d). $$
%
% $$ $$ Note that since Matlab starts indexing at $1$, we need to shift $r$ by one.
% $$ $$
%%
clear all % this clears all variables of their values

%% Numerical rank
% Before we start, one should recall the trouble one can have with rank
% estimates. In the second line, the value 1e-16 might be the result of
% round-off errors or simply correct. 
rank(diag([1,1e-15]))
rank(diag([1,1e-16]))

%% TT-SVD
% The following procedure constructs a TT-representation for any input
% tensor T. Since the input tensor may have low rank, we require a
% truncation tolerace _tol_ , in order to avoid the trouble indicated above
% (this is exactly what the rank function does, but we are going to do this ourselves). 
%
% There is no need yet to fully understand this computation. In the next 
% exercise, we only want to whether it is correct. For more insight, have
% a look at the worksheet _Hierarchy in the game of modes_.

tol = 1e-15;
T1 = randn(4,3,4,2);

n = size(T1)
d = length(n)
G = cell(1,d);
r = zeros(1,d+1);
r(1) = 1; r(d+1) = 1; % note the shift

B = T1;
for mu = 1:d
    [U,S,V] = svd(reshape(B,[n(mu)*r(mu),prod(n(mu+1:d))]),0);
    r(mu+1) = sum(diag(S) > tol); % count singular values larger than tol
    G{mu} = permute( reshape(U,[r(mu),n(mu),r(mu+1)]) , [1,3,2]);
    B = S*V';
end
G{d} = G{d} * B; % B at this point is a scalar
r
%%
size(G{2})
[r(2),r(3),n(2)]
G{2}(:,:,2) % this corresponds to G_2(2)

%% EXERCISE 1: entrywise evaluation
% Complete the following function eval_TT_entry. As input, it expects
% a TT-representation G of a tensor T and an array I. The output is the
% value of the entry T_I. As before, it is more convenient to
% hand over n as parameter as well.
%
I = [1,2,3,2]
T1(1,2,3,2)
eval_TT_entry(G,I,n)

%% Core product
% We can treat the cores as individual objects and multiply them:

%%
%   function H = boxtimes(G1,G2)
% 
%   [r1,k1,n1] = size(G1);
%   [r2,k2,n2] = size(G2);
% 
%   H = zeros(r1,k2,n1,n2);
%   for i = 1:n1
%       for j = 1:n2
%           H(:,:,i,j) = G1(:,:,i) * G2(:,:,j);
%       end
%   end
%   H = reshape(H,[r1,k2,n1*n2]);
% 
%   end
H = boxtimes(G{2},G{3});
[r(2),r(4),n(2)*n(3)]
size(H)
%%
H_unfold = reshape(H,[r(2),r(4),n(2),n(3)]);
H_unfold(:,:,2,3)
G{2}(:,:,2) * G{3}(:,:,3)
H(:,:,8)

%%
% We can also build the full tensor:
T1_test = boxtimes(boxtimes(boxtimes(G{1},G{2}),G{3}),G{4});
size(T1_test)
T1_test = reshape(T1_test,n);
T1_test(:,:,2,2)
T1(:,:,2,2)

%% EXERCISE 2: About the core product:
% Explain the outputs above. Why is _boxtimes_ defined this way?

%% Left and right interface matrcies
% Another important aspect is a certain reshaping of matrices.

%%
%   function R = right_fold(Gmu)
%       [r1,r2,n] = size(Gmu);
%       R = reshape(Gmu,[r1,r2*n]);
%   end
%

%%
%   function L = left_fold(Gmu) 
%       [r1,r2,n] = size(Gmu);
%       L = reshape(  permute(Gmu,[1,3,2])  ,[r1*n,r2]);
%   end
Gmu = zeros(2,3,2); 
Gmu(:) = 1:2*3*2
left_fold(Gmu)
right_fold(Gmu)

%%
% The inverse operations are as follows:
Gmu
%%
right_unfold(right_fold(Gmu),2,3,2)
%%
left_unfold(left_fold(Gmu),2,3,2)

%% Left- and right-orthogonality
% Like in the Tucker format, orthogonality constraints play an important
% role. Since the TT format is ordered, one core is either
% left or right of another one, like when multiplying matrices.
% Analogously, cores can either be left- or right-orthogonal, which
% corresponds to column and row orthogonality. Due to the contruction
% above, G{1} to G{d-1} are left-orthogonal, but not G{d}.
round10 = @(x) round(x*1e10)/1e10;
L = left_fold(G{1});
round10(L'*L)
%%
L = left_fold(G{2});
round10(L'*L)
%%
L = left_fold(G{3});
round10(L'*L)
%%
L = left_fold(G{4});
round10(L'*L)

%% 
% For the norm of the tensor then follows
norm(T1(:))
norm(G{4}(:))

%%
% Sometimes one wishes a specific distribution of orthogonality contraints:
G34 = boxtimes(G{3},G{4});

R = right_fold(G{4});
[Q,R] = qr(R',0); Q = Q'; R = R';

G{4} = right_unfold(Q,r(4),r(5),n(4));
G{3} = boxtimes(G{3},R);

G34_new = boxtimes(G{3},G{4});
round10(norm(G34(:)-G34_new(:))) % we have not changed the product of these cores

%%
% Hence now...
L = left_fold(G{3});
round10(L'*L)
%%
R = right_fold(G{3});
round10(R*R')
%%
R = right_fold(G{4});
round10(R*R')

%%
% If we now want to know the norm of T1, we have to use G{3}
norm(G{3}(:))
norm(G{4}(:))
norm(T1(:)) % Frobenius norm

%% EXERCISE 3: norm of a TT represented tensor
% Given a randomly initialized representation, calculate the Frobenius norm
% of the tensor being represented, without constructing the full tensor.
%
% To validate your result, you may then construct the full tensor.
d = 4
G = cell(1,d);
n = randi([3,4],1,d)
r = [1,randi([2,3],1,d-1),1]
for mu = 1:d
    G{mu} = randn(r(mu),r(mu+1),n(mu));
end

%% SOLUTION 3
'ANSWER 3 MISSING';

%% EXERCISE 4: comparision of CP, Tucker and TT format
% Recall the properties of the CP and Tucker format and compare them to the
% TT format. What is a weakspot of TT, which both CP and Tucker do not
% have? What about tensors with very large dimension? 

%% EXERCISE* 5: compression within a TT representation
% Create a function in a separat Matlab file which applies the following
% procedure to a tensor, given by a TT-representation, 
% but without constructing the full tensor:
%
%%
%   function [T,r] = TT_truncate_full_tensor(T,tol)
% 
%   n = size(T);
%   d = length(n);
%   r = ones(1,d+1);
% 
%   for mu = 1:d-1
%       [U,S,V] = svd( reshape(T,[prod(n(1:mu)),prod(n(mu+1:d))]) ,0);
%       r(mu+1) = sum(diag(S)>tol); % count number of entries larger than tol
%       T = U(:,1:(mu+1)) * S(1:(mu+1),1:(mu+1)) * V(:,1:(mu+1))';
%   end
%   
%   T = reshape(T,n);
%   end
%
d = randi([5,6],1)
G = cell(1,d);
n = randi([4,5],1,d)
r = [1,randi([2,4],1,d-1),1]
%%
for mu = 1:d % random representation
    G{mu} = randn(r(mu),r(mu+1),n(mu));
    G{mu} = boxtimes(G{mu},diag(5.^(-(1:r(mu+1)))));
end

T2 = 1; % construct full tensor
for mu = 1:d
    T2 = boxtimes(T2,G{mu});
end
T2 = reshape(T2,n);
%%
[T2_trunc,r] = TT_truncate_full_tensor(T2,1e-4);
T2_trunc(1:16)
r

%% SOLUTION 5
'ANSWER 5 MISSING';





























