function [t,est] = lra(A, m, varargin)
[matvec,n] = process_matrix(A, varargin{:});

addpath('../code')
S = generate_test_matrix(n,floor(m/2),'signs',varargin{:});
[Q,~] = qr(matvec(S),0);
t = trace(Q'*matvec(Q));
est = NaN;

end