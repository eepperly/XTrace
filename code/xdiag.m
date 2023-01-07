function d = xdiag(matvec, adjvec, n, m, varargin)
m = floor(m/2);

%% Choose test matrix
[Om,improved] = generate_test_matrix(n,m,'signs',varargin{:});
if improved
    error('Improved test vectors not implemented for XDiag')
end

%% Randomized SVD
Y = matvec(Om); 
[Q,R] = qr(Y,0);

%% Quantities needed for diagonal estimation
Z = adjvec(Q); T=Z'*Om;
S = cnormr(eye(size(R,1)) - triu(R,1) / R)';
dQZ = diag_prod(Q',Z'); dQSSZ = diag_prod((Q*S)',(Z*S)');
dOmQT = diag_prod(Om',(Q*T).'); dOmY = diag_prod(Om',Y.');
dOmQSST = diag_prod(Om',(Q*S*diag(diag_prod(S,T))).');

%% Diagonal estimate
d = dQZ + (-dQSSZ+dOmY-dOmQT+dOmQSST)/m;
end