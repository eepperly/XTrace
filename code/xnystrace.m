function [t,err] = xnystrace(matvec, n, m, varargin)
%% Choose test matrix
[Om,improved] = generate_test_matrix(n,m,'signs',varargin{:});

%% Nystrom approximation
Y = matvec(Om);
nu = eps*norm(Y,'fro')/sqrt(n);
Y = Y + nu*Om; % Shift for numerical stability
[Q,R] = qr(Y,0);
H = Om'*Y; C = chol((H+H')/2);
B = R/C; % Nystrom approx is Q*B*B'*Q'

%% Normalization
if improved
    [QQ,RR] = qr(Om,0);
    WW = QQ'*Om; 
    SS = cnormr(eye(size(RR,1)) - triu(RR,1) / RR)';
    scale = (n - m + 1) ./ (n - vecnorm(WW).^2 ...
        + (abs(diag_prod(SS,WW)') .* vecnorm(SS)).^2);
else
    scale = ones(1,m);
end

% Trace estimate
warning('off','MATLAB:nearlySingularMatrix');
W = Q'*Om; S = (B/C') .* (diag(inv(H))').^(-0.5);
warning('on','MATLAB:nearlySingularMatrix');
dSW = diag_prod(S, W).'; dOmY = diag_prod(Om, Y).';
ests = norm(B,'fro')^2 - vecnorm(S).^2 + (dOmY - vecnorm(B'*W).^2 ...
    + abs(dSW).^2) .* scale - nu*n;
t = mean(ests);
err = std(ests)/sqrt(m);
end
