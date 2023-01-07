function A = rand_with_evals(evals,varargin)
if ~isempty(varargin)
    complex = varargin{1};
else
    complex = 0;
end
n = length(evals);
[Q,~] = qr(randn(n) + 1i*randn(n)*complex);
A = Q * diag(evals) * Q';
end

