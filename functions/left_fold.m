function L = left_fold(Gmu)
[r1,r2,n] = size(Gmu);

L = reshape(permute(...
    Gmu,[1,3,2]),...
    [r1*n,r2]);
end