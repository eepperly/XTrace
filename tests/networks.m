%% Setup
addpath('../code')
addpath('../existing_estimators')
check_expmv
rng(42)

if exist('yeast.mat', 'file')
    load('yeast.mat')
else
    error(['To run this script, must have yeast.mat in the current' ...
        'folder. Can be downloaded at: ' ...
        'https://sparse.tamu.edu/Pajek/yeast'])
end
A = Problem.A;
d = eig(A);
trials = 1000;
close all

%% Methods
methods = {@bks, @diag_plusplus, @xdiag};
method_names = {'Hutch', 'Diag++', 'XDiag'};
colors = {[0 0.4470 0.7410], "#EDB120", [0.4940 0.1840 0.5560]};
markers = 'o*x';

ms = 10:10:200;

for i = 1:2
    if i == 1
        target = diag(expm(A));
        matvec = @(x) expmv(1,A,x);
        algs = [1 2 3];
    elseif i == 2
        target = diag(A^3)/2;
        matvec = @(x) A*(A*(A*x)) / 2;
        algs = [1 3];
    end

    for j = algs
        method = methods{j};
        
        errors = [];
        for m = ms
            m
            errs = [];
            ests = [];
            for l = 1:trials
                d = method(matvec, m, size(A,1), 'signs');
                errs(end+1) = norm(d - target,"inf");
            end
            errors(:,end+1) = mean(errs) / norm(target,"inf");
        end
        
        figure(i)
        semilogy(ms, errors, 'LineWidth', 1.5, 'Color',...
            colors{j}, 'Marker', markers(j), ...
            'MarkerSize', 10, 'MarkerFaceColor', colors{j});
        hold on
        set(gca, 'YScale', 'log')
        drawnow
    end

    set(gca, 'FontSize', 20)
    if i == 1
        legend({'BKS','Diag++','XDiag'},'FontSize',20)
        saveas(gcf,'../figs/ee_diag.png');
        saveas(gcf,'../figs/ee_diag.fig');
    elseif i == 2
        saveas(gcf,'../figs/triangle_diag.png');
        saveas(gcf,'../figs/triangle_diag.fig');
    end
end