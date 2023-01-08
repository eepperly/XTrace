%% Setup
addpath('../code')
addpath('../existing_estimators')
rng(42)

%% Generate test matrices

n = 1000;
As = {rand_with_evals(linspace(1,3,n)),...
      rand_with_evals((1:n).^(-2)),...
      rand_with_evals(0.7.^(0:(n-1))),...
      rand_with_evals([ones(ceil(n/20),1);1e-3*ones(n-ceil(n/20),1)])};
names = {'flat','poly','fastexp','smallstep'}; 
ranges = {[1e-2 1e0], [1e-5 1e0], [1e-10 1e0], [2e-4 1e0]};

%% Parameters for experiments

num_trials = 1;
num_range = 20;

%% Run experiments
close all
for i = 1:length(As)
    A = As{i};
    trace_A = trace(A);
    name = names{i};
    range = ranges{i};
    range = logspace(log10(range(1)),log10(range(2)),num_range);
    
    ahpp_errs = zeros(num_range,1);
    ms = zeros(num_range,1);
    errs_run = zeros(num_trials,1);
    ms_run = zeros(num_trials,1);
    for k = 1:num_range
        fprintf('%s A-Hutch++ %e\n', name, range(k))
        for l = 1:num_trials
            if mod(l,floor(num_trials/10)) == 0
                fprintf('.')
            end
            [t,ms_run(l)] = adap_hpp(n,@(x) A*x,range(k)*trace_A,0.1);
            errs_run(l) = abs(trace_A - t) / trace_A; 
        end
        ahpp_errs(k) = mean(errs_run);
        ms(k) = mean(ms_run);
        ahpp_errs(k)
        ms(k)
        fprintf('\n')
    end

    ms_xtrace = round(linspace(min(ms),min(max(ms),1000),num_range));
    errors_xtrace = zeros(size(ms_xtrace));
    for k = 1:length(ms_xtrace)
        m = ms_xtrace(k);
        fprintf('%s XTrace %d\n', name, m)
        for l = 1:num_trials
            if mod(l,floor(num_trials/10)) == 0
                fprintf('.')
            end
            errors_xtrace(k) = errors_xtrace(k) ...
                + abs(trace_A - xtrace(@(x) A*x,n,m,'signs'))...
                /trace_A/num_trials;
        end
        fprintf('\n')
    end
    
    figure(i)
    semilogy(ms, ahpp_errs, 'LineWidth', 1.5, 'MarkerSize', 10, ...
            'Color', "#A2142F", 'MarkerFaceColor', "#A2142F",...
            'Marker', 'diamond');
    hold on
    semilogy(ms, range, 'k--', 'LineWidth', 1.5);
    semilogy(ms_xtrace, errors_xtrace, 'LineWidth', 1.5,...
        'MarkerSize', 10, 'Color', "#7E2F8E", 'MarkerFaceColor',...
        "#7E2F8E", 'Marker', 'x');
    xlabel('$m$','FontSize',24)
    ylabel('Average relative error','FontSize',24)
    drawnow
    saveas(gcf, sprintf('../figs/%s_adapt.png', name))
    saveas(gcf, sprintf('../figs/%s_adapt.fig', name))
end