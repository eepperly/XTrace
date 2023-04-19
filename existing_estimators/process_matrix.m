function [matvec,n,adjvec,m] = process_matrix(A,varargin)
if isnumeric(A)
    matvec = @(x) A*x;
    n = size(A,1);
    adjvec = @(x) A'*x;
    m = size(A,2);
    return
end

if ~isa(A, 'function_handle')
    error('Input must be either a numeric matrix or a function handle');
end
matvec = A;

n_set = false;
m_set = false;
adjvec_set = false;
for i = 1:length(varargin)
    if ~n_set && isscalar(varargin{i})
        n = varargin{i}; n_set = true;
        m = n;
    elseif ~m_set && isscalar(varargin{i})
        m = varargin{i}; m_set = true;
    elseif ~adjvec_set && isa(varargin{i}, 'function_handle')
        adjvec = varargin{i}; adjvec_set = true;
    end

    if n_set && adjvec_set
        break
    end
end

if ~adjvec_set; adjvec = matvec; end
if ~n_set
    error(['Size of input matrix must be specified if matrix'...
        ' is provided as a function handle'])
end
end