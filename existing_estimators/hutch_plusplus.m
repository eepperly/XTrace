function [t,est] = hutch_plusplus(matvec, n, m, varargin)
% HUTCH_PLUSPLUS Hutch++ trace estimator
% Implementation is a modification of:
%    https://github.com/RaphaelArkadyMeyerNYU/HutchPlusPlus/blob/main/simple/simple_hutchplusplus.m

addpath('../code')
S = generate_test_matrix(n,ceil(m/3),'signs',varargin{:});
G = generate_test_matrix(n,floor(m/3),'signs',varargin{:});

[Q,~] = qr(matvec(S),0);
G = G - Q*(Q'*G);

t = trace(Q'*matvec(Q)) + 1/size(G,2)*trace(G'*matvec(G));
est = NaN;
end