function [t,est] = hutch(matvec, n, m, varargin)

addpath('../code')
Omega = generate_test_matrix(n,m,'signs',varargin{:});
t = trace(Omega'*matvec(Omega)) / m;
est = NaN;

end

