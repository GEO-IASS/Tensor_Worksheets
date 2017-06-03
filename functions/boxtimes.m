function H = boxtimes(G1,G2)

[r1,k1,n1] = size(G1);
[r2,k2,n2] = size(G2);

H = zeros(r1,k2,n1,n2);
for i = 1:n1
    for j = 1:n2
        H(:,:,i,j) = G1(:,:,i) * G2(:,:,j);
    end
end
H = reshape(H,[r1,k2,n1*n2]);

end