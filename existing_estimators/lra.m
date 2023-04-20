function [t,est,time,processtime] = lra(A, m, varargin)
[matvec,n] = process_matrix(A, varargin{:});
tic;
S = generate_test_matrix(n,floor(m/2),'signs',varargin{:});
processtime = toc;
tic; Y = matvec(S); matvectime = toc;
tic; [Q,~] = qr(Y,0); processtime = processtime + toc;
tic; Y = matvec(Q); matvectime = matvectime + toc;
tic; t = trace(Q'*Y); processtime = processtime + toc;
est = NaN;
time = processtime + matvectime;
end