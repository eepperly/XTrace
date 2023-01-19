function [t,est] = hutch(A, m, varargin)
[matvec,n] = process_matrix(A, varargin{:});
Omega = generate_test_matrix(n,m,'signs',varargin{:});
t = trace(Omega'*matvec(Omega)) / m;
est = NaN;
end

