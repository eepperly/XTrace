function [Om,improved] = generate_test_matrix(n,m,default,varargin)
if isempty(varargin)
    type = default;
else
    type = varargin{1};
end

improved = false;
if strcmp(type, 'improved')
    improved = true;
    Om = sqrt(n) * cnormc(randn(n,m));
elseif strcmp(type, 'cimproved')
    improved = true;
    Om = sqrt(n) * cnormc(randn(n,m) + 1i*randn(n,m));
elseif strcmp(type, 'rademacher') || strcmp(type, 'signs')
    Om = -3 + 2*randi(2,n,m);
elseif strcmp(type, 'steinhaus') || strcmp(type, 'csigns') ...
    || strcmp(type, 'phases')
    Om = exp(2*pi*1i*rand(n,m)); 
elseif strcmp(type, 'gaussian')
    Om = randn(n,m);
elseif strcmp(type, 'cgaussian')
    Om = 1/sqrt(2) * randn(n,m) + 1i/sqrt(2)*randn(n,m);
elseif strcmp(type, 'unif') || strcmp(type, 'sphere')
    Om = sqrt(n) * normc(randn(n,m));
elseif strcmp(type, 'cunif') || strcmp(type, 'csphere')
    Om = sqrt(n) * normc(randn(n,m)+1i*randn(n,m));
elseif strcmp(type, 'orth')
    Om = sqrt(n)*orth(randn(n,m));
elseif strcmp(type, 'corth')
    Om = sqrt(n)*orth(randn(n,m) + 1i*randn(n,m));
else 
    error('"%s" not recognized as matrix type', type)
end
end