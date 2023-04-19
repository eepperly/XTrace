function [t,err,time,processtime] = xnystrace(A, m, varargin)
[matvec,n] = process_matrix(A, varargin{:});

%% Generate test matrix
tic;
[Om,improved] = generate_test_matrix(n,m,'improved',varargin{:});
processtime = toc;

%% Matrix-vector products
tic;
Y = matvec(Om);
matvectime = toc;

%% Finish computation
tic;
[t,err] = xnystrace_helper(Y, Om, improved);
processtime = processtime + toc;
time = processtime + matvectime;
end

