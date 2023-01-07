function M = cnormc(M)
M = M ./ vecnorm(M,2,1);
end

