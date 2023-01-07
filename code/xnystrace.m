function [t,err] = xnystrace(matvec, n, m, varargin)
%% Choose test matrix
[Om,improved] = generate_test_matrix(n,m,'signs',varargin{:});
if improved
    error('Improved test vectors not implemented for XNysTrace')
end

%% Nystrom approximation
Y = matvec(Om);
nu = eps*norm(Y,'fro')/sqrt(n);
Y = Y + nu*Om; % Shift for numerical stability
[Q,R] = qr(Y,0);
H = Om'*Y; C = chol((H+H')/2);
B = R/C; % Nystrom approx is Q*B*B'*Q'

%% Normalization
if improved
    W = Q'*Om; 
    S = cnormr(eye(size(R,1)) - triu(R,1) / R)';
    scale = sqrt(n - m + 1) ./ sqrt(n - vecnorm(W).^2 ...
        + (abs(diag_prod(S,W)') .* vecnorm(S)).^2);
    Om = Om .* scale;
end

% Trace estimate
warning('off','MATLAB:nearlySingularMatrix');
W = Q'*Om; S = (B/C') .* (diag(inv(H))').^(-0.5);
warning('on','MATLAB:nearlySingularMatrix');
dSW = diag_prod(S, W).'; dOmY = diag_prod(Om, Y).';
ests = norm(B,'fro')^2 - vecnorm(S).^2 + dOmY - vecnorm(B'*W).^2 ...
    + abs(dSW).^2 - nu*n;
t = mean(ests);
err = std(ests)/sqrt(m);
end

