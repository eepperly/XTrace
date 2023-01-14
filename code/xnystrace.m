function [t,err] = xnystrace(matvec, n, m, varargin)
[Om,improved] = generate_test_matrix(n,m,'improved',varargin{:});
Y = matvec(Om);
[t,err] = xnystrace_helper(Y, Om, improved);
end

