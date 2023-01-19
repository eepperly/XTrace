addpath('../code')
check_expmv
rng(42)

hs = logspace(-1,1,10);
betas = logspace(-1,1,10);
n = 18;
tol = 1e-4;

Zs = zeros(length(hs),length(betas));
Z_ms = zeros(length(hs),length(betas));
Z_ests = zeros(length(hs),length(betas));
Es = zeros(length(hs),length(betas));
EZ_ms = zeros(length(hs),length(betas));
EZ_ests = zeros(length(hs),length(betas));

for i = 1:length(hs)
    h = hs(i);
    H = tfim(n,h);
    for j = 1:length(betas)
        fprintf('%d / %d\n',(i-1)*length(hs) + j,length(hs)*length(betas));
        beta = betas(j);
        b = -(1+h)*n;
        H_shift = H - b*speye(2^n,2^n);
        fprintf('**Z**\n')
        [Zs(i,j),Z_ests(i,j),Z_ms(i,j)] = ...
            xnystrace_tol(@(x) expmv(1,-beta*H_shift,x),...
            0,tol,size(H,1),'signs');
        fprintf('**EZ**\n')
        [EZ,EZ_ests(i,j),EZ_ms(i,j)] = ...
            xnystrace_tol(@(x) expmv(1,-beta*H_shift,H_shift*x), ...
            0,tol,size(H,1),'signs');
        Es(i,j) = EZ / Zs(i,j) + b;
    end
end

[xx,yy] = meshgrid(betas,hs);
contourf(xx,yy,Es/n,-logspace(-1,1.2,12)); set(gca,'ColorScale','log');
colorbar('Ticks',[-10,-3,-1,-0.3,-0.1])
set(gca,'XScale', 'log')
set(gca,'YScale', 'log')