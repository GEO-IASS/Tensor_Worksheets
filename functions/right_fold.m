function R = right_fold(Gmu)
[r1,r2,n] = size(Gmu);

R = reshape(Gmu,[r1,r2*n]);
end