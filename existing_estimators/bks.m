function [d,err,time,processtime] = bks(A,m,varargin)
[matvec,n] = process_matrix(A, varargin{:});
tic; Om = generate_test_matrix(n,m,'signs',varargin{:}); processtime = toc;
tic; Y = matvec(Om); matvectime = toc;
tic; d = diag_prod(Om',Y.') / m;
processtime = processtime + toc; time = processtime + matvectime;
err = NaN(size(d));
end