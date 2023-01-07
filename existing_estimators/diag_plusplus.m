function d = diag_plusplus(matvec, adjvec, n, m, varargin)
% HUTCH_PLUSPLUS Hutch++ trace estimator
% Implementation is a modification of:
%    https://github.com/RaphaelArkadyMeyerNYU/HutchPlusPlus/blob/main/simple/simple_hutchplusplus.m

addpath('../code')
S = generate_test_matrix(n,ceil(m/3),'signs',varargin{:});
G = generate_test_matrix(n,floor(m/3),'signs',varargin{:});

[Q,~] = qr(matvec(S),0);
B = matvec(G - Q*(Q'*G));

d = diag_prod(Q',(Q'*matvec(Q))*Q') ...
    + sum(G .* (B - Q*(Q'*B))) / size(G,2);
    
end