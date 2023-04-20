function [t,est,time,processtime] = hutch(A, m, varargin)
[matvec,n] = process_matrix(A, varargin{:});
tic;
Omega = generate_test_matrix(n,m,'signs',varargin{:});
processtime = toc;
tic; Y = matvec(Omega); matvectime = toc;
tic; t = trace(Omega'*Y) / m; processtime = processtime + toc;
est = NaN;
time = processtime + matvectime;
end

