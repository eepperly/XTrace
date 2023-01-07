function M = cnormr(M)
M = M ./ vecnorm(M,2,2);
end

