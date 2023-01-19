%% Set up

rng(42)
addpath('../code/')
addpath('../existing_estimators/')

%% Generate test matrices

n = 1000;
As = {rand_with_evals(linspace(1,3,n)),...
      rand_with_evals((1:n).^(-2)),...
      rand_with_evals(0.9.^(0:(n-1))),...
      rand_with_evals(0.7.^(0:(n-1))),...
      rand_with_evals([ones(ceil(n/20),1);1e-3*ones(n-ceil(n/20),1)]),...
      rand_with_evals([ones(ceil(n/20),1);1e-8*ones(n-ceil(n/20),1)])};
names = {'flat','poly','slowexp','fastexp','smallstep','bigstep'}; 

%% Methods

methods = {@hutch, @lra, @hutch_plusplus, @nystrompp, @xtrace, @xnystrace};
method_names = {'Hutch', 'LRA', 'Hutch++', 'Nystr\"om++', 'XTrace',...
    'NysTrace'};
colors = {"#0072BD","#D95319","#EDB120","#4DBEEE","#7E2F8E","#77AC30"};
markers = 'os*<x^';

%% Parameters for experiments

ms = 20:20:300;
num_trials = 1000;

%% Run experiments
close all
for i = 1:length(As)
    A = As{i};
    trace_A = trace(A);
    name = names{i};
    
    for j = 1:length(methods)
        method = methods{j};
        color = colors{j};
        marker = markers(j);
        
        errors = [];
        for m = ms
            errs = [];
            for l = 1:num_trials
                errs(end+1) = abs(method(A, m, 'signs') - trace_A);
            end
            errors(end+1) = mean(errs) / trace_A;
        end
        
        figure(i)
        semilogy(ms, errors, 'LineWidth', 1.5, 'MarkerSize', 10, ...
            'Color', color, 'MarkerFaceColor', color, 'Marker', marker);
        hold on
        set(gca, 'YScale', 'log')
        drawnow
    end
    xlabel('$m$','FontSize',24)
    ylabel('Average relative error','FontSize',24)
    if i == 1
        legend(method_names,'Location','Best')
    end
    saveas(gcf, sprintf('../figs/%s.png', name))
    saveas(gcf, sprintf('../figs/%s.fig', name))
end