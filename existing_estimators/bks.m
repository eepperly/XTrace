function d = bks(matvec,adjvec,n,m,varargin)
addpath('../code')
Omega = generate_test_matrix(n,m,'signs',varargin{:});
d = diag_prod(Om',matvec(Om).') / m;
end