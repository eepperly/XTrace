function [t,err] = xtrace(A, m, varargin)
[matvec,n] = process_matrix(A, varargin{:});
m = floor(m/2);
[Om,improved] = generate_test_matrix(n,m,'improved',varargin{:});
Y = matvec(Om);
[Q,R] = qr(Y,0);
Z = matvec(Q);
[t,err] = xtrace_helper(Om,Z,Q,R,improved);
end