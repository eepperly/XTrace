function d = diag_prod(A,B)
d = sum(conj(A).*B,1).'; % computes diag(A'*B)
end