function [t,est] = lra(matvec, n, m, varargin)

addpath('../code')
S = generate_test_matrix(n,floor(m/2),'signs',varargin{:});
[Q,~] = qr(matvec(S),0);
t = trace(Q'*matvec(Q));
est = NaN;

end