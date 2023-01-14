function [t,err] = xnystrace_helper(Y, Om, improved)
%% Nystrom 
[n,m] = size(Y);
nu = eps*norm(Y,'fro')/sqrt(n);
Y = Y + nu*Om; % Shift for numerical stability
[Q,R] = qr(Y,0);
H = Om'*Y; C = chol((H+H')/2);
B = R/C; % Nystrom approx is Q*B*B'*Q'

%% Normalization
if improved
    [QQ,RR] = qr(Om,0);
    WW = QQ'*Om; 
    warning('off','MATLAB:nearlySingularMatrix');
    SS = cnormc(inv(RR)');
    warning('on','MATLAB:nearlySingularMatrix');
    scale = (n - m + 1) ./ (n - vecnorm(WW).^2 ...
        + abs(diag_prod(SS,WW)' .* vecnorm(SS)).^2);
else
    scale = ones(1,m);
end

% Trace estimate
warning('off','MATLAB:nearlySingularMatrix');
W = Q'*Om; S = (B/C') .* (diag(inv(H))').^(-0.5);
warning('on','MATLAB:nearlySingularMatrix');
dSW = diag_prod(S, W).';
ests = norm(B,'fro')^2 - vecnorm(S).^2 + abs(dSW).^2 .* scale - nu*n;
t = mean(ests);
err = std(ests)/sqrt(m);
end

