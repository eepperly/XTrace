function [t,err,m] = xtrace_tol(A, abstol, reltol, varargin)
[matvec,n] = process_matrix(A, varargin{:});
Y = zeros(n,0); Om = zeros(n,0); Z = zeros(n,0);
err = Inf;
if (reltol ~=0); t = Inf; else; t = 0; end
while err >= abstol + reltol * abs(t)
    % Randomized SVD    
    [NewOm,improved] = generate_test_matrix(n,max(2,size(Y,2)), ...
        'improved',varargin{:});
    Y = [Y matvec(NewOm)]; Om = [Om NewOm];%#ok<AGROW> 
    [Q,R] = qr(Y,0);
    Z = [Z matvec(Q(:,size(Z,2)+1:end))]; %#ok<AGROW> 
    [t,err] = xtrace_helper(Om, Z, Q, R, improved);
end
m = 2*size(Z,2);
end