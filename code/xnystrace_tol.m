function [t,err,m] = xnystrace_tol(matvec, n, abstol, reltol, varargin)
Y = zeros(n,0); Om = zeros(n,0);
err = Inf;
if (reltol ~=0); t = Inf; else; t = 0; end
while err >= abstol + reltol * t
    % Nystrom approximation
    [NewOm,improved] = generate_test_matrix(n,max(2,size(Y,2)), ...
        'improved',varargin{:});
    Y = [Y matvec(NewOm)]; Om = [Om NewOm];%#ok<AGROW>
    [t,err] = xnystrace_helper(Y, Om, improved);
end
m = size(Y,2);
end

