function trest = nystrompp(A,m,varargin)
% NYSTROMPP Nystrom++ trace estimator
% Implementation is a modification of:
%   https://github.com/davpersson/A-Hutch-/blob/main/Nystrom%2B%2B/nystrompp.m
% which carries the following license:
% BSD 2-Clause License
% 
% Copyright (c) 2021, David Persson
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
[matvec,n] = process_matrix(A, varargin{:});

%Generate random matrices
Omega = generate_test_matrix(n,ceil(m/2),'signs',varargin{:});
Psi = generate_test_matrix(n,floor(m/2),'signs',varargin{:});

%Compute matrix products with A
Y = matvec(Omega); Z = matvec(Psi);

%Regularizaton
nu = sqrt(n)*eps(norm(Y));
Y_nu = Y + nu*Omega;
C = chol(Omega'*Y_nu);
B = Y_nu/C;
[U,S,~] = svd(B,'econ');
Lambda = max(0,S^2-nu*eye(size(S)));

tr1 = trace(Lambda);
tr2 = 2*(product_trace(Psi',Z) - product_trace(Psi',U*(Lambda*(U'*Psi))))/m;
trest = tr1 + tr2;

end

function t = product_trace(A,B)
%Computes trace(A*B) without computing all entries of A*B
t = sum(sum(A.*B',2));
end