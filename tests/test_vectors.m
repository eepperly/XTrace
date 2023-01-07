%% Set up
addpath('../code')
rng(42)

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
test_vector_types = {'rademacher', 'gaussian', 'unif', 'improved'};
markers = {'x','v','>','pentagram'};
markersizes = [7 10 15 10];
styles = {':','-','--','-.'};
colors = {"#0072BD","#D95319","#EDB120","#7E2F8E"};
markercolors = {"#0072BD","#D95319",'none',"#7E2F8E"};

%% Parameters for experiments

ks = 20:20:300;
num_trials = 10;

%% Run experiments
close all
for i = 1:length(As)
    A = As{i};
    trace_A = trace(A);
    name = names{i};
    
    for j = 1:length(test_vector_types)
        test_vector_type = test_vector_types{j};
        
        errors = [];
        stds = [];
        for k = ks
            errs = [];
            for l = 1:num_trials
                errs(end+1) = abs(xtrace(@(x) A*x, n, k, ...
                    test_vector_type) - trace_A);
            end
            errors(end+1) = mean(errs) / trace_A;
            stds(end+1) = std(errs) / trace_A;
        end
        
        figure(i)
        semilogy(ks, errors, 'LineWidth', 1.5, 'Marker', markers{j},...
            'Color', colors{j}, 'MarkerFaceColor', markercolors{j},...
            'LineStyle', styles{j}, 'MarkerSize', markersizes(j)); hold on
        set(gca, 'YScale', 'log')
        drawnow
    end
    xlabel('Matrix--vector products $m$','FontSize',24)
    ylabel('Average relative error','FontSize',24)
    if i == 1
        legend(test_vector_types,'Location','Best')
    end
    saveas(gcf, sprintf('../figs/test_vector_%s.png', name))
    saveas(gcf, sprintf('../figs/test_vector_%s.fig', name))
end