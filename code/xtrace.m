function [t,err,time,processtime] = xtrace(A, m, varargin)
[matvec,n] = process_matrix(A, varargin{:});
m = floor(m/2);

%% Generate test matrix
tic;
[Om,improved] = generate_test_matrix(n,m,'improved',varargin{:});
processtime = toc;

%% First round of matvecs
tic;
Y = matvec(Om);
matvectime = toc;

%% Orthonormalize
tic;
[Q,R] = qr(Y,0);
processtime = processtime + toc;

%% First round of matvecs
tic;
Z = matvec(Q);
matvectime = matvectime + toc;

%% Complete
tic;
[t,err] = xtrace_helper(Om,Z,Q,R,improved);
processtime = processtime + toc;
time = processtime + matvectime;
end