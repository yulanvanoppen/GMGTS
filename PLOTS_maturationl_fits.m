clearvars

files = dir(fullfile('experimental', '*.mat'));
files = {files.name};


for index = 1:length(files)
    
load(['experimental/' files{index}])


%%
close all

figure('position', [70, 10, 800, 600])
tiledlayout(2, 3);

legends = [];
parameter_names = system.parameters_variable;

labels = {'Data', 'GMGTS', 'GTS'};

CI_GMGTS = quantile(estimator.results_GMGTS.fitted_ss, [.05 .95], 3);
% CI_GTS = quantile(estimator.results_GTS.fitted2, [.05 .95], 3);

nexttile(1)
title("Random Effects")
hold on
beta_GMGTS = estimator.results_GMGTS.beta_fs;
% beta_GTS = estimator.results_GTS.beta_fs;
b_GMGTS = estimator.results_GMGTS.b_est;
% b_GTS = estimator.results_GTS.b_est;
D_GMGTS = estimator.results_GMGTS.D_est;
% D_GTS = estimator.results_GTS.D_GTS;
            
h3 = scatter(beta_GMGTS(:, 1), beta_GMGTS(:, 2), [], [0, 0.4470, 0.7410], 'square', 'MarkerEdgeAlpha', .3, 'LineWidth', 1);
% h4 = scatter(beta_GTS(:, 1), beta_GTS(:, 2), [], [0.4940, 0.1840, 0.5560], 'square', 'MarkerEdgeAlpha', .3, 'LineWidth', 1);

h1 = Estimator.mvlncontour(b_GMGTS, D_GMGTS, '^', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1);
% h2 = Estimator.mvncontour(b_GTS, D_GTS, 'square', 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 1);

xlabel('$\mathsf k_\textsf{p}$', 'interpreter', 'latex')
ylabel('$\mathsf k_\textsf{m}$', 'interpreter', 'latex')
xlim([max(0, exp(b_GMGTS(1)-3*sqrt(D_GMGTS(1, 1)))) exp(b_GMGTS(1)+3*sqrt(D_GMGTS(1, 1)))])
ylim([max(0, exp(b_GMGTS(2)-3*sqrt(D_GMGTS(2, 2)))) exp(b_GMGTS(2)+3*sqrt(D_GMGTS(2, 2)))])
ylim(ylim + [0 -.1*range(ylim)])
legend([h3 h1], {'Estimates', 'Mean', 'Covariance'}, 'Location', 'northeast');



nexttile(4)
title("Dark Protein $\boldsymbol D$", 'interpreter', 'latex')
xlabel('Time (min)')
ylabel('Measurement (A.U.)')
hold on
            
h1 = plot(estimator.results_GMGTS.t_fine, estimator.results_GMGTS.population(:, 1), ...
          ':', 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 1, 1), ':', 'LineWidth', 1.1, 'Color', [0, 0.4470, 0.7410]);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 1, 2), ':', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
h2 = patch([estimator.results_GMGTS.t_fine flip(estimator.results_GMGTS.t_fine)], [CI_GMGTS(:, 1, 1)' flip(CI_GMGTS(:, 1, 2))'], ...
           [0, 0.4470, 0.7410], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');

if system.K == 3
    h3 = plot(estimator.results_GMGTS.t_fine, estimator.results_GMGTS.population(:, 2), ...
              ':', 'LineWidth', 2, 'Color', [0.3010 0.7450 0.9330]);
    plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 2, 1), ':', 'LineWidth', 1.1, 'Color', [0.3010 0.7450 0.9330]);
    plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 2, 2), ':', 'LineWidth', 1.5, 'Color', [0.3010 0.7450 0.9330]);
    h4 = patch([estimator.results_GMGTS.t_fine flip(estimator.results_GMGTS.t_fine)], [CI_GMGTS(:, 2, 1)' flip(CI_GMGTS(:, 2, 2))'], ...
               [0.3010 0.7450 0.9330], 'EdgeAlpha', 0, 'FaceAlpha', .4, 'LineStyle', 'none');
           
    legend([h2 h4], {'$D_1$', '$D_2$'}, 'Location', 'northwest', 'interpreter', 'latex')
end
            
% h1 = plot(estimator.results_GTS.t, estimator.results_GTS.population(:, 1), ...
%           '-.', 'LineWidth', 2, 'Color', [0.4940, 0.1840, 0.5560]);
% plot(estimator.results_GTS.t, CI_GTS(:, 1, 1), '-.', 'LineWidth', 1.5, 'Color', [0.4940, 0.1840, 0.5560]);
% plot(estimator.results_GTS.t, CI_GTS(:, 1, 2), '-.', 'LineWidth', 1.5, 'Color', [0.4940, 0.1840, 0.5560]);
% h2 = patch([estimator.results_GTS.t flip(estimator.results_GTS.t)], [CI_GTS(:, 1, 1)' flip(CI_GTS(:, 1, 2))'], ...
%            [0.4940, 0.1840, 0.5560], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');

ylim(max(0, ylim))


nexttile([2 2])
title("Fluorescent Protein $\boldsymbol F$", 'interpreter', 'latex')
xlabel('Time (min)')
ylabel('Measurement (A.U.)')
hold on

h = plot(estimator.results_GMGTS.t_data, ...
    reshape(estimator.results_GMGTS.traces(:, 1, 1:estimator.results_GMGTS.N), [], estimator.results_GMGTS.N), '-');
set(h, {'color'}, num2cell(parula(estimator.results_GMGTS.N), 2));
            
h0 = plot(estimator.results_GMGTS.t_fine, estimator.results_GMGTS.population(:, end), ...
          '--', 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, end, 1), ':', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, end, 2), ':', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
h1 = patch([estimator.results_GMGTS.t_fine flip(estimator.results_GMGTS.t_fine)], [CI_GMGTS(:, end, 1)' flip(CI_GMGTS(:, end, 2))'], ...
           [0, 0.4470, 0.7410], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');
            
% h0 = plot(estimator.results_GTS.t, estimator.results_GTS.population(:, end), ...
%           '-.', 'LineWidth', 2, 'Color', [0.4940, 0.1840, 0.5560]);
% plot(estimator.results_GTS.t, CI_GTS(:, end, 1), '-.', 'LineWidth', 1.5, 'Color', [0.4940, 0.1840, 0.5560]);
% plot(estimator.results_GTS.t, CI_GTS(:, end, 2), '-.', 'LineWidth', 1.5, 'Color', [0.4940, 0.1840, 0.5560]);
% h2 = patch([estimator.results_GTS.t flip(estimator.results_GTS.t)], [CI_GTS(:, end, 1)' flip(CI_GTS(:, end, 2))'], ...
%            [0.4940, 0.1840, 0.5560], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');

ylim(max(0, ylim))
legend([h(1) h0 h1], {'Measurements', 'Predicted process mean', 'Predicted process distribution'}, 'Location', 'northwest')


%%
filestr = sprintf('figures/maturation_%s.tex', data(idx).name);
matlab2tikz(filestr)

fid = fopen(filestr, 'r');
f = fread(fid,'*char')';
f = strrep(f, 'only marks, ', '');
f = strrep(f, 'mark=square, ', 'only marks, mark=square, opacity=.5, ');
f = strrep(f, 'mark=triangle, ', 'only marks, mark=triangle, ');
f = strrep(f, 'xlabel style={', 'xlabel style={at={(axis description cs: 0.5, 0.01)}, ');
f = strrep(f, 'ylabel style={', 'ylabel style={at={(axis description cs: 0.03, 0.5)}, ');
f = strrep(f, 'contour/draw color=mycolor1', 'contour/draw color=mycolor1, color=mycolor1');
f = strrep(f, 'mark size=1.0607pt', 'mark size=2pt');
fclose(fid);
fid = fopen(filestr, 'w');
fwrite(fid, f);
fclose(fid);

end











