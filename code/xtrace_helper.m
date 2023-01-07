function [t,err] = xtrace_helper(Y, Om, Z, Q, R, improved)

%% Normalization
[n,m] = size(Om);
W = Q'*Om; 
S = cnormr(eye(size(R,1)) - triu(R,1) / R)';
if improved
    scale = (n - m + 1) ./ (n - (vecnorm(W)').^2 ...
        + (abs(diag_prod(S,W)) .* (vecnorm(S)')).^2);
else
    scale = ones(m,1);
end

%% Quantities needed for trace estimation
H = Q'*Z; HW = H*W; T=Z'*Om;
dSW = diag_prod(W, S); dSHS = diag_prod(S, H*S);
dOmY = diag_prod(Om, Y); dWR = diag_prod(W, R);
dTW = diag_prod(T, W); dWHW = diag_prod(W, HW);
dSRmHW = diag_prod(S, R-HW); dTmHRS = diag_prod(T-H'*W,S);

% Trace estimate
ests = trace(H)*ones(m,1) - dSHS + (dOmY - dWR - dTW + dWHW...
    + conj(dSW) .* dSRmHW + abs(dSW).^2 .* dSHS + dTmHRS .* dSW) .* scale;
t = mean(ests);
err = std(ests)/sqrt(m);
end