function [d,err,time,processtime] = diag_plusplus(A, m, varargin)
% DIAG_PLUSPLUS Diag++ trace estimator
% Implementation is a modification of:
%    https://github.com/RaphaelArkadyMeyerNYU/HutchPlusPlus/blob/main/simple/simple_hutchplusplus.m
[matvec,n] = process_matrix(A, varargin{:});

tic;
S = generate_test_matrix(n,ceil(m/3),'signs',varargin{:});
G = generate_test_matrix(n,floor(m/3),'signs',varargin{:});
processtime = toc; 

tic; Y = matvec(S); matvectime = toc;
tic
[Q,~] = qr(Y,0);
Y = G - Q*(Q'*G);
processtime = processtime + toc;
tic;
B = matvec(Y);
matvectime = matvectime + toc;

tic;
d = diag_prod(Q',(Q'*matvec(Q))*Q') ...
    + sum(G .* (B - Q*(Q'*B))) / size(G,2);
processtime = processtime + toc;
err = NaN(size(d));
time = processtime + matvectime;
end