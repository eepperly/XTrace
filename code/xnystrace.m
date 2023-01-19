function [t,err] = xnystrace(A, m, varargin)
[matvec,n] = process_matrix(A, varargin{:});
[Om,improved] = generate_test_matrix(n,m,'improved',varargin{:});
Y = matvec(Om);
[t,err] = xnystrace_helper(Y, Om, improved);
end

