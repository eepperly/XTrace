function d = xdiag(A, m, varargin)
[matvec,n,adjvec] = process_matrix(A, varargin{:});
m = floor(m/2);

%% Choose test matrix
[Om,~,type] = generate_test_matrix(n,m,'signs',varargin{:});
if ~strcmp(type, 'signs')
    error('XDiag only implemented with random signs')
end

%% Randomized SVD
Y = matvec(Om); 
[Q,R] = qr(Y,0);

%% Quantities needed for diagonal estimation
Z = adjvec(Q); T=Z'*Om;
warning('off','MATLAB:nearlySingularMatrix');
S = cnormc(inv(R)');
warning('on','MATLAB:nearlySingularMatrix');
dQZ = diag_prod(Q',Z'); dQSSZ = diag_prod((Q*S)',(Z*S)');
dOmQT = diag_prod(Om',(Q*T).'); dOmY = diag_prod(Om',Y.');
dOmQSST = diag_prod(Om',(Q*S*diag(diag_prod(S,T))).');

%% Diagonal estimate
d = dQZ + (-dQSSZ+dOmY-dOmQT+dOmQSST)/m;
end