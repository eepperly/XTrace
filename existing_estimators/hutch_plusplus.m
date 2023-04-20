function [t,est,time,processtime] = hutch_plusplus(A, m, varargin)
% HUTCH_PLUSPLUS Hutch++ trace estimator
% Implementation is a modification of:
%    https://github.com/RaphaelArkadyMeyerNYU/HutchPlusPlus/blob/main/simple/simple_hutchplusplus.m
[matvec,n] = process_matrix(A, varargin{:});

tic;
S = generate_test_matrix(n,ceil(m/3),'signs',varargin{:});
G = generate_test_matrix(n,floor(m/3),'signs',varargin{:});
processtime = toc;

tic; Y = matvec(S); matvectime = toc;
tic;
[Q,~] = qr(Y,0);
G = G - Q*(Q'*G);
processtime = processtime + toc;

tic; Y = matvec(Q); Z = matvec(G); matvectime = matvectime + toc;
tic;
t = trace(Q'*Y) + 1/size(G,2)*trace(G'*Z);
processtime = processtime + toc;
est = NaN;
time = matvectime + processtime;
end