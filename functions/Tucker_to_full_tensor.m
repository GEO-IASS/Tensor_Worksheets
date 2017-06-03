function C = Tucker_to_full_tensor(d,C,U)

for mu = 1:d
    C = mu_mode_prod(C,d,U{mu},mu);
end

end
