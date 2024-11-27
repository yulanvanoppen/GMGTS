clearvars

% files = dir(fullfile('simulation', '*.mat'));
% files = {files.name};

load(['simulation/maturation'])


%%
close all

figure('position', [70, 10, 800, 600])
tl = tiledlayout(system.P, system.P);

legends = [];
parameter_names = system.parameters_variable;

labels = {'GMGTS.I', 'GMGTS.II', 'GTS.I', 'GTS.II'};
labels2 = {'GMGTS', 'GMGTS', 'Truth'};

for p = 1:system.P
    for q = 1:p
        nexttile(q + (p-1)*system.P)
        hold on
        if q == p
            beta_GMGTS = estimator.results_GMGTS.beta_fs;
            beta_GTS = estimator.results_GTS.beta_fs;
            beta_true = ground_truth.beta;
            b_GMGTS = estimator.results_GMGTS.b_est;
            b_GTS = estimator.results_GTS.b_est;
            b_true = ground_truth.b;
            D_GMGTS = estimator.results_GMGTS.D_GTS;
            D_GTS = estimator.results_GTS.D_GTS;
            D_true = ground_truth.D;
            
            h1 = histogram(beta_GMGTS(:, p), 20, 'Normalization', 'pdf', ...
                           'FaceColor', [0, 0.4470, 0.7410], 'FaceAlpha', .3);
                       
            grid = linspace(min(beta_GMGTS(:, p)), max(beta_GMGTS(:, p)), 100);
            pdf = normpdf(grid, b_GMGTS(p), sqrt(D_GMGTS(p, p)));
            h2 = plot(grid, pdf, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5);
            patch([grid flip(grid)], [pdf 0*pdf], [0, 0.4470, 0.7410], 'EdgeAlpha', 0, 'FaceAlpha', .4);
            
            h3 = histogram(beta_GTS(:, p), 20, 'Normalization', 'pdf', ...
                           'FaceColor', [0.4940, 0.1840, 0.5560], 'FaceAlpha', .3);
                       
            grid = linspace(min(beta_GTS(:, p)) - range(beta_GTS(:, p))/4, max(beta_GTS(:, p)) + range(beta_GTS(:, p))/4, 100);
            pdf = normpdf(grid, b_GTS(p), sqrt(D_GTS(p, p)));
            h4 = plot(grid, pdf, 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 1.5);
            patch([grid flip(grid)], [pdf 0*pdf], [0.4940, 0.1840, 0.5560], 'EdgeAlpha', 0, 'FaceAlpha', .4);
            
            h5 = histogram(beta_true(:, p), 20, 'Normalization', 'pdf', ...
                           'FaceColor', [0.8500, 0.3250, 0.0980], 'FaceAlpha', .3);
                       
            grid = linspace(min(beta_true(:, p)) - range(beta_true(:, p))/4, max(beta_true(:, p)) + range(beta_true(:, p))/4, 100);
            pdf = normpdf(grid, b_true(p), sqrt(D_true(p, p)));
            h6 = plot(grid, pdf, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5);
            patch([grid flip(grid)], [pdf 0*pdf], [0.8500, 0.3250, 0.0980], 'EdgeAlpha', 0, 'FaceAlpha', .4);
            
            if p == 2
                legend([h1 h3 h5], labels2);
            end
            
        else
            h1 = Estimator.mvncontour(b_GMGTS([q p]), D_GMGTS([q p], [q p]), '^', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1);
            h2 = Estimator.mvncontour(b_GTS([q p]), D_GTS([q p], [q p]), 'square', 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 1);
            h3 = Estimator.mvncontour(b_true([q p]), D_true([q p], [q p]), 'pentagram', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1);
            
            xl = xlim;
            yl = ylim;
            
            h4 = scatter(beta_GMGTS(:, q), beta_GMGTS(:, p), [], [0, 0.4470, 0.7410], '^', 'MarkerEdgeAlpha', .3, 'LineWidth', 1);
            h5 = scatter(beta_GTS(:, q), beta_GTS(:, p), [], [0.4940, 0.1840, 0.5560], 'square', 'MarkerEdgeAlpha', .3, 'LineWidth', 1);
            h6 = scatter(beta_true(:, q), beta_true(:, p), [], [0.8500, 0.3250, 0.0980], 'pentagram', 'MarkerEdgeAlpha', .3, 'LineWidth', 1);
            
            xlim([max(0, xl(1)) xl(2)])
            ylim([max(0, yl(1)) yl(2)])
%             if p == system.P && q == 2
%                 legend([h1(2) h2(2) h3(2)], labels2, 'Location', 'northeast');
%             end
        end
        if q == 1, xlim([.01 .1]), end
        if q == 2, xlim([5e-3 1e-1]), end
        if q == 3, xlim([0 .75]), end
        if p == system.P
            xlabel(parameter_names(q))
            if q ~= 3, ylim([0, .3]), end
        end
    end
end

xlim([0, .3])
% switch marker
%     case '^', col = [0, 0.4470, 0.7410]; center = ':v'; lw = 1.2;
%     case 'pentagram', col = [0.8500, 0.3250, 0.0980]; center = '--hexagram'; lw = 1;
%     case 'square', col = [0.4940, 0.1840, 0.5560]; center = '-.diamond'; lw = 1;
% end

 

%%
filestr = '../../writing/manuscript/figures/maturation_simulation_distribution.tex';
matlab2tikz(filestr)

% fid = fopen(filestr, 'r');
% f = fread(fid,'*char')';
% f = strrep(f, 'xminorticks=true,', ...
%            'xminorticks=true, minor xtick={.5, .6, .7, .8, .9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120},');
% f = strrep(f, 'xmode=log,', 'xmode=log, xticklabels={{.5}, {1}, {2}, {5}, {10}, {20}, {50}, {100}},');
% fclose(fid);
% fid = fopen(filestr, 'w');
% fwrite(fid, f);
% fclose(fid);




