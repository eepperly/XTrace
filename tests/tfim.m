function H = tfim(n,h)
N = 2^n;
rows = zeros(N*(n+1), 1); cols = zeros(N*(n+1), 1);
vals = zeros(N*(n+1), 1); idx = 1;

for i = uint64(0):uint64(N-1)
    for j = 0:(n-1)
        mask = uint64(2^j);
        rows(idx) = i+1; cols(idx) = bitxor(i,mask)+1;
        vals(idx) = -h; idx = idx+1;
    end
    val = bitxor(bitget(i,1),bitget(i,n));
    for j = 1:(n-1)
        val = val + bitxor(bitget(i,j),bitget(i,j+1));
    end
    rows(idx) = i+1; cols(idx) = i+1;
    vals(idx) = 2*double(val)-n; idx = idx+1;
end
H = sparse(rows, cols, vals, N, N);
end