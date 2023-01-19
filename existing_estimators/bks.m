function d = bks(A,m,varargin)
[matvec,n] = process_matrix(A, varargin{:});
addpath('../code')
Om = generate_test_matrix(n,m,'signs',varargin{:});
d = diag_prod(Om',matvec(Om).') / m;
end