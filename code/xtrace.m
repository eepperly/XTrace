function [t,err] = xtrace(matvec, n, m, varargin)
m = floor(m/2);

%% Choose test matrix
[Om,improved] = generate_test_matrix(n,m,'improved',varargin{:});

%% Randomized SVD
Y = matvec(Om); 
[Q,R] = qr(Y,0);

%% Normalization
W = Q'*Om; 
S = cnormr(eye(size(R,1)) - triu(R,1) / R)';
if improved
    scale = sqrt(n - m + 1) ./ sqrt(n - vecnorm(W).^2 ...
        + (abs(diag_prod(S,W)') .* vecnorm(S)).^2);
    Om = Om .* scale; R = R .* scale; Y = Y .* scale; W = W .* scale;
end

%% Quantities needed for trace estimation
Z = matvec(Q); H = Q'*Z; HW = H*W; T=Z'*Om;
dSW = diag_prod(W, S); dSHS = diag_prod(S, H*S);
dOmY = diag_prod(Om, Y); dWR = diag_prod(W, R);
dTW = diag_prod(T, W); dWHW = diag_prod(W, HW);
dSRmHW = diag_prod(S, R-HW); dTmHRS = diag_prod(T-H'*W,S);

% Trace estimate
ests = trace(H)*ones(m,1) - dSHS + dOmY - dWR - dTW + dWHW...
    + conj(dSW) .* dSRmHW + abs(dSW).^2 .* dSHS + dTmHRS .* dSW;
t = mean(ests);
err = std(ests)/sqrt(m);
end