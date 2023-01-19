function [matvec,n,adjvec] = process_matrix(A,varargin)
if isnumeric(A)
    matvec = @(x) A*x;
    n = size(A,1);
    adjvec = @(x) A'*x;
    return
end

if ~isa(A, 'function_handle')
    error('Input must be either a numeric matrix or a function handle');
end
matvec = A;

n_set = false;
adjvec_set = false;
for i = 1:length(varargin)
    if ~n_set && isscalar(varargin{i})
        n = varargin{i}; n_set = true;
    elseif ~adjvec_set && isa(varargin{i}, 'function_handle')
        adjvec = varargin{i}; adjvec_set = true;
    end

    if n_set && adjvec_set
        break
    end
end

if ~adjvec_set; adjvec = matvec; end
end