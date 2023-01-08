function [trest,total_matvecs,lowrank_matvecs,trest_matvecs] = adap_hpp(matrix_size,Afun,epsilon,delta)
% This is a direct copy of 
%   https://github.com/davpersson/A-Hutch-/blob/main/A-Hutch%2B%2B/adap_hpp.m
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

%INPUTS
%A: Square matrix
%delta: Failure probability. Algorithm succeeds wpral 1-2*delta
%epsilon: Error tolerance

%OUTPUTS:
%trest: Approximation to trace(A)
%matvecs: Total number of matvecs used
%rsvd: Total number of matvecs used for randomized svd
%M: Total number of matvecs used for Hutchinson estimation

%Set up parameters
C = (4*log(2/delta)/epsilon^2); %Constant for number of matvecs

flag = 1;
fnc = [];

%First iteration
y = Afun(randn(matrix_size,1));
q = y/norm(y);
Q = q;

x = Afun(q);
t = q'*x;
c = t^2;
trest1 = t;
iteration = 1;
b = norm(x)^2;
fnc(iteration) = 2*iteration+C*(c - 2*b);

%Remaining iterations
while flag
    
    %Get new column in Q
    y = Afun(randn(matrix_size,1));
    qt = y-Q*(Q'*y); 
    
    if norm(qt) < 1e-10
        
        lowrank_matvecs = 2*iteration;
        trest_matvecs = 0;
        total_matvecs = lowrank_matvecs + trest_matvecs;
        trest = trest1;
        return
        
        
    end
    
    qt = qt/norm(qt); qt = qt - Q*(Q'*qt); q = qt/norm(qt);
    Q = [Q q];
    x = Afun(q);
    
    %Update recursion
    b = b + norm(x)^2;
    t = q'*x;
    trest1 = trest1 + t;
    c = c + 2*norm(Q(:,1:iteration)'*x)^2 + t^2;
    
    %Update function
    iteration = iteration + 1;
    fnc(iteration) = 2*iteration+C*(c - 2*b);
    
    %Check if iteration should be stopped
    if (iteration > 2) && (fnc(iteration-1) < fnc(iteration)) && (fnc(iteration-2) < fnc(iteration-1))
        
        flag = 0;
        
    end
    
end

lowrank_matvecs = 2*iteration;

%Combine Hutchinson and Frobenius norm estimation
iteration = 0;
flag = 1;
t = 0;
trest2 = [];

while flag
    
    psi = randn(matrix_size,1);
    y = psi - Q*(Q'*psi);
    y = Afun(y);
    y = y - Q*(Q'*y);
    trest2 = [trest2 psi'*y];
    t = t + y'*y;
    estFrob = t/(iteration+1);
    alpha = supfind(iteration+1,delta);
    M = ceil(C*estFrob/alpha);
    
    if iteration + 1 > M
        
        flag = 0;
        
    end
    
    iteration = iteration + 1;
    
end

%Perform Hutchinson estimation
trest2 = mean(trest2);
trest_matvecs = iteration;


%Return outputs
trest = trest1 + trest2;
total_matvecs = lowrank_matvecs + trest_matvecs; 
end