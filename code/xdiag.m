function [d,est,time,processtime] = xdiag(A, m, varargin)
[matvec,n,adjvec] = process_matrix(A, varargin{:});
m = floor(m/2);

%% Choose test matrix
tic;
[Om,~,type] = generate_test_matrix(n,m,'signs',varargin{:});
if ~strcmp(type, 'signs')
    error('XDiag only implemented with random signs')
end
processtime = toc;

%% Randomized SVD
tic; Y = matvec(Om); matvectime = toc;
tic; [Q,R] = qr(Y,0); processtime = processtime + toc;

%% Quantities needed for diagonal estimation
tic; Z = adjvec(Q); matvectime = matvectime + toc;
tic;
T=Z'*Om;
warning('off','MATLAB:nearlySingularMatrix');
S = cnormc(inv(R)');
warning('on','MATLAB:nearlySingularMatrix');
dQZ = diag_prod(Q',Z'); dQSSZ = diag_prod((Q*S)',(Z*S)');
dOmQT = diag_prod(Om',(Q*T).'); dOmY = diag_prod(Om',Y.');
dOmQSST = diag_prod(Om',(Q*S*diag(diag_prod(S,T))).');

%% Diagonal estimate
d = dQZ + (-dQSSZ+dOmY-dOmQT+dOmQSST)/m;
processtime = processtime + toc;
time = processtime + matvectime;
est = NaN(size(d));
end