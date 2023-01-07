function d = tfim_eigs(n,h)
if mod(n,2) ~= 0
    error('Not implemented for n odd')
end
ks = ((-n+1):2:(n-1))*pi/n;
energies = 2*sqrt(1+h^2+2*h*cos(ks)).';
ground = -0.5*sum(energies);
mat = bitand(uint64(0:(2^n-1))',uint64(2.^(0:(n-1)))) ~= 0;
d = mat(mod(sum(mat,2),2)==0,:)*energies + ground;

ks = ((-n/2):(n/2-1))*2*pi/n;
energies = 2*sqrt(1+h^2+2*h*cos(ks)).';
energies(1) = -2*(1+h);
energies(n/2+1) = 2*(1-h);
ground = -0.5*sum(energies);
d = [d;mat(mod(sum(mat,2),2)==1,:)*energies + ground];
end