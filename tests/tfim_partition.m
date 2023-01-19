%% Add code
addpath('../code')
addpath('../existing_estimators')
check_expmv
rng(42)

%% Run experiments

n = 18; h = 1; beta = 3;
H = tfim(n,h);
d = tfim_eigs(n,h);
b = -(1+h)*n;
H_shift = H - b*speye(2^n,2^n);
exact = sum(exp(-beta*(d-b)));
ms = round(10 * 2.^(0:(1/3):4));
num_trials = 100;

hutch_errs = zeros(length(ms),num_trials);
xnystrace_errs = zeros(length(ms),num_trials);
xnystrace_ests = zeros(length(ms),num_trials);
xtrace_errs = zeros(length(ms),num_trials);
xtrace_ests = zeros(length(ms),num_trials);
hutchpp_errs = zeros(length(ms),num_trials);
for i = 1:length(ms)
    m = ms(i);
    fprintf('%d\t',m)
    for j = 1:num_trials
        fprintf('.')
        t = hutch(@(x) expmv(1,-beta*H_shift,x),m,2^n);
        hutch_errs(i,j) = abs(t - exact) / exact;
        [t,s] = xnystrace(@(x) expmv(1,-beta*H_shift,x),m,2^n,'signs');
        xnystrace_ests(i,j) = s / exact;
        xnystrace_errs(i,j) = abs(t - exact) / exact;
        [t,s] = xtrace(@(x) expmv(1,-beta*H_shift,x),m,2^n,'signs');
        xtrace_ests(i,j) = s / exact;
        xtrace_errs(i,j) = abs(t - exact) / exact;
        t = hutch_plusplus(@(x) expmv(1,-beta*H_shift,x),m,2^n,'signs');
        hutchpp_errs(i,j) = abs(t - exact) / exact;
    end
    fprintf('\n')
end

%% Plots
figure(1)
loglog(ms,mean(hutch_errs,2),'Color',"#0072BD",'Marker','o', ...
                'MarkerSize', 10, 'MarkerFaceColor', '#0072BD'); hold on
loglog(ms,mean(hutchpp_errs,2),'Color',"#EDB120",'Marker','*', ...
                'MarkerSize', 10);
loglog(ms,mean(xtrace_errs,2),'Color','#7E2F8E','Marker','x', ...
                'MarkerSize', 10)
loglog(ms,mean(xnystrace_errs,2),'Color','#77AC30','Marker','^', ...
                'MarkerSize', 10, 'MarkerFaceColor', '#77AC30')
axis([1e1 Inf 1e-15 1e0])
xlabel('Matrix--vector products $m$')
ylabel('Mean relative error')
set(gca,'FontSize',20)
legend({'Hutch','Hutch++','XTrace','XNysTrace'},'FontSize',20,...
    'Location','best')
saveas(gcf,'../figs/tfim_errs_all.png')
saveas(gcf,'../figs/tfim_errs_all.fig')

figure(2)
loglog(ms,mean(xtrace_errs,2),'Color','#7E2F8E','Marker','x', ...
                'MarkerSize', 10); hold on
loglog(ms,mean(xtrace_ests,2),'--','Color','#7E2F8E')
loglog(ms,mean(xnystrace_errs,2),'Color','#77AC30','Marker','^', ...
                'MarkerSize', 10, 'MarkerFaceColor', '#77AC30')
loglog(ms,mean(xnystrace_ests,2),'--','Color','#77AC30')
axis([1e1 Inf 1e-15 1e-5])
xlabel('Matrix--vector products $m$')
ylabel('Mean relative error')
set(gca,'FontSize',20)
legend({'XTrace','','XNysTrace',''},'FontSize',20,...
    'Location','best')
saveas(gcf,'../figs/tfim_errs_with_ests.png')
saveas(gcf,'../figs/tfim_errs_with_ests.fig')