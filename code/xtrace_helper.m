function [t,err] = xtrace_helper(Om, Z, Q, R, improved)

%% Normalization
[n,m] = size(Om);
W = Q'*Om; 
warning('off','MATLAB:nearlySingularMatrix');
S = cnormc(inv(R)');
warning('on','MATLAB:nearlySingularMatrix');
if improved
    scale = (n - m + 1) ./ (n - (vecnorm(W)').^2 ...
        + abs(diag_prod(S,W) .* vecnorm(S)').^2);
else
    scale = ones(m,1);
end

%% Quantities needed for trace estimation
H = Q'*Z; HW = H*W; T=Z'*Om;
dSW = diag_prod(S, W); dSHS = diag_prod(S, H*S);
dTW = diag_prod(T, W); dWHW = diag_prod(W, HW);
dSRmHW = diag_prod(S, R-HW); dTmHRS = diag_prod(T-H'*W,S);

% Trace estimate
ests = trace(H)*ones(m,1) - dSHS + (- dTW + dWHW...
    + conj(dSW) .* dSRmHW + abs(dSW).^2 .* dSHS + dTmHRS .* dSW) .* scale;
t = mean(ests);
err = std(ests)/sqrt(m);
end