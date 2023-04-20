%% Add code
addpath('../code')
addpath('../existing_estimators')
check_expmv
rng(42)

%% Run experiments

n = 18; h = 1; beta = 0.6;
H = tfim(n,h);
d = tfim_eigs(n,h);
b = -(1+h)*n;
H_shift = H - b*speye(2^n,2^n);
exact = sum(exp(-beta*(d-b)));
ms = round(10 * 2.^(0:(1/3):4));
num_trials = 10;

hutch_errs = zeros(length(ms),num_trials);
hutch_times = zeros(length(ms),num_trials);
hutch_preps = zeros(length(ms),num_trials);
xnystrace_errs = zeros(length(ms),num_trials);
xnystrace_times = zeros(length(ms),num_trials);
xnystrace_preps = zeros(length(ms),num_trials);
xtrace_errs = zeros(length(ms),num_trials);
xtrace_times = zeros(length(ms),num_trials);
xtrace_preps = zeros(length(ms),num_trials);
hutchpp_errs = zeros(length(ms),num_trials);
hutchpp_times = zeros(length(ms),num_trials);
hutchpp_preps = zeros(length(ms),num_trials);
for i = 1:length(ms)
    m = ms(i);
    fprintf('%d\t',m)
    for j = 1:num_trials
        fprintf('.')
        [t,~,hutch_times(i,j),hutch_preps(i,j)] ...
            = hutch(@(x) expmv(1,-beta*H_shift,x),m,2^n,'signs');
        hutch_errs(i,j) = abs(t - exact) / exact;
        [t,~,xnystrace_times(i,j),xnystrace_preps(i,j)] ...
            = xnystrace(@(x) expmv(1,-beta*H_shift,x),m,2^n,'signs');
        xnystrace_errs(i,j) = abs(t - exact) / exact;
        [t,~,xtrace_times(i,j),xtrace_preps(i,j)] ...
            = xtrace(@(x) expmv(1,-beta*H_shift,x),m,2^n,'signs');
        xtrace_errs(i,j) = abs(t - exact) / exact;
        [t,~,hutchpp_times(i,j),hutchpp_preps(i,j)]...
            = hutch_plusplus(@(x) expmv(1,-beta*H_shift,x),m,2^n,'signs');
        hutchpp_errs(i,j) = abs(t - exact) / exact;
    end
    fprintf('\n')
end

%% Plots
figure(1)
loglog(mean(hutch_preps,2),mean(hutch_errs,2),'Color',"#0072BD",'Marker','o', ...
                'MarkerSize', 10, 'MarkerFaceColor', '#0072BD'); hold on
loglog(mean(hutchpp_preps,2),mean(hutchpp_errs,2),'Color',"#EDB120",'Marker','*', ...
                'MarkerSize', 10);
loglog(mean(xtrace_preps,2),mean(xtrace_errs,2),'Color','#7E2F8E','Marker','x', ...
                'MarkerSize', 10)
loglog(mean(xnystrace_preps,2),mean(xnystrace_errs,2),'Color','#77AC30','Marker','^', ...
                'MarkerSize', 10, 'MarkerFaceColor', '#77AC30')
axis([1e1 Inf 1e-15 1e0])
xlabel('Processing time (sec)')
ylabel('Mean relative error')
set(gca,'FontSize',20)
legend({'Hutch','Hutch++','XTrace','XNysTrace'},'FontSize',20,...
    'Location','best')
saveas(gcf,'../figs/process_timing.png')
saveas(gcf,'../figs/process_timing.fig')

figure(1)
loglog(mean(hutch_times,2),mean(hutch_errs,2),'Color',"#0072BD",'Marker','o', ...
                'MarkerSize', 10, 'MarkerFaceColor', '#0072BD'); hold on
loglog(mean(hutchpp_times,2),mean(hutchpp_errs,2),'Color',"#EDB120",'Marker','*', ...
                'MarkerSize', 10);
loglog(mean(xtrace_times,2),mean(xtrace_errs,2),'Color','#7E2F8E','Marker','x', ...
                'MarkerSize', 10)
loglog(mean(xnystrace_times,2),mean(xnystrace_errs,2),'Color','#77AC30','Marker','^', ...
                'MarkerSize', 10, 'MarkerFaceColor', '#77AC30')
axis([1e1 Inf 1e-15 1e0])
xlabel('Total time (matvecs + processing, sec)')
ylabel('Mean relative error')
set(gca,'FontSize',20)
saveas(gcf,'../figs/total_timing.png')
saveas(gcf,'../figs/total_timing.fig')