clearvars
    
load('simulation/maturation_accuracy2.mat')


%%
close all

figure('position', [70, 10, 800, 600])
tiledlayout(2, 3);

legends = [];
parameter_names = system.parameters_variable;

labels = {'Data', 'GMGTS', 'GTS'};

CI_GMGTS = quantile(estimator.results_GMGTS.fitted2, [.05 .95], 3);
CI_GTS = quantile(estimator.results_GTS.fitted2, [.05 .95], 3);
CI_true = quantile(ground_truth.original, [.05 .95], 3);

colGMGTS = [0, 0.4470, 0.7410];
colGTS   = [0.9290 0.6940 0.1250];
coltrue  = [0.8500 0.3250 0.0980];

nexttile(1)
title("Random Effects")
hold on
beta_GMGTS = estimator.results_GMGTS.beta_fs;
beta_GTS = estimator.results_GTS.beta_fs;
beta_true = ground_truth.beta;
b_GMGTS = estimator.results_GMGTS.b_est;
b_GTS = estimator.results_GTS.b_est;
b_true = ground_truth.b;
D_GMGTS = estimator.results_GMGTS.D_GTS;
D_GTS = estimator.results_GTS.D_GTS;
D_true = ground_truth.D;
            
% h3 = scatter(beta_GMGTS(:, 1), beta_GMGTS(:, 2), [], colGMGTS, 'square', 'MarkerEdgeAlpha', .3, 'LineWidth', 1);
% h4 = scatter(beta_GTS(:, 1), beta_GTS(:, 2), [], colGTS, 'square', 'MarkerEdgeAlpha', .3, 'LineWidth', 1);

h1 = Estimator.mvncontour(b_true, D_true, 'o', 'Color', coltrue, 'LineWidth', 1);
h2 = Estimator.mvncontour(b_GTS, D_GTS, 'square', 'Color', colGTS, 'LineWidth', 1);
h3 = Estimator.mvncontour(b_GMGTS, D_GMGTS, '^', 'Color', colGMGTS, 'LineWidth', 1);

xlabel('$k_\t{p}$', 'interpreter', 'latex')
ylabel('$k_\t{m}$', 'interpreter', 'latex')
ylim([max(0, min(b_GMGTS(1)-3*sqrt(D_GMGTS(1, 1)), b_GTS(2)-3*sqrt(D_GTS(1, 1)))) ...
      max(b_GMGTS(1)+3*sqrt(D_GMGTS(1, 1)), b_GTS(1)+3*sqrt(D_GTS(1, 1)))])
ylim([max(0, min(b_GMGTS(2)-3*sqrt(D_GMGTS(2, 2)), b_GTS(2)-3*sqrt(D_GTS(2, 2)))) ...
      max(b_GMGTS(2)+3*sqrt(D_GMGTS(2, 2)), b_GTS(2)+3*sqrt(D_GTS(2, 2)))])
ylim(ylim + [0 .5*range(ylim)])
legend([h1(2) h2(2) h3(2)], {'Truth', 'GTS', 'GMGTS'}, 'Location', 'northeast');



nexttile(4)
title("Dark Protein $\boldsymbol{D}$", 'interpreter', 'latex')
xlabel('Time (min)')
ylabel('Concentration (A.U.)')
hold on
            
h1 = plot(ground_truth.t, ground_truth.moriginal(:, 1), ...
          '-.', 'LineWidth', 2, 'Color', colGTS);
plot(ground_truth.t, CI_true(:, 1, 1), '-.', 'LineWidth', 1.5, 'Color', coltrue);
plot(ground_truth.t, CI_true(:, 1, 2), '-.', 'LineWidth', 1.5, 'Color', coltrue);
h2 = patch([ground_truth.t flip(estimator.results_GTS.t)], [CI_true(:, 1, 1)' flip(CI_true(:, 1, 2))'], ...
           coltrue, 'EdgeAlpha', 0, 'FaceAlpha', .15, 'LineStyle', 'none');
       
h1 = plot(estimator.results_GTS.t, estimator.results_GTS.population(:, 1), ...
          '-.', 'LineWidth', 2, 'Color', colGTS);
plot(estimator.results_GTS.t, CI_GTS(:, 1, 1), '-.', 'LineWidth', 1.5, 'Color', colGTS);
plot(estimator.results_GTS.t, CI_GTS(:, 1, 2), '-.', 'LineWidth', 1.5, 'Color', colGTS);
h2 = patch([estimator.results_GTS.t flip(estimator.results_GTS.t)], [CI_GTS(:, 1, 1)' flip(CI_GTS(:, 1, 2))'], ...
           colGTS, 'EdgeAlpha', 0, 'FaceAlpha', .15, 'LineStyle', 'none');
       
h1 = plot(estimator.results_GMGTS.t_fine, estimator.results_GMGTS.population(:, 1), ...
          ':', 'LineWidth', 2, 'Color', colGMGTS);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 1, 1), ':', 'LineWidth', 1.1, 'Color', colGMGTS);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 1, 2), ':', 'LineWidth', 1.5, 'Color', colGMGTS);
h2 = patch([estimator.results_GMGTS.t_fine flip(estimator.results_GMGTS.t_fine)], [CI_GMGTS(:, 1, 1)' flip(CI_GMGTS(:, 1, 2))'], ...
           colGMGTS, 'EdgeAlpha', 0, 'FaceAlpha', .15, 'LineStyle', 'none');
            
ylim(max(0, ylim))




nexttile([2 2])
title("Fluorescent Protein $\boldsymbol{F}$", 'interpreter', 'latex')
xlabel('Time (min)')
ylabel('Concentration (A.U.)')
hold on

h = plot(estimator.results_GMGTS.t_data, ...
    reshape(estimator.results_GMGTS.traces(:, 1, 1:estimator.results_GMGTS.N), [], estimator.results_GMGTS.N), '-', 'LineWidth', .9);
set(h, {'color'}, num2cell(repmat(mean(parula(estimator.results_GMGTS.N), 2), 1, 3), 2));
            
h1 = plot(ground_truth.t, ground_truth.moriginal(:, 2), ...
          '-.', 'LineWidth', 2, 'Color', colGTS);
plot(ground_truth.t, CI_true(:, 2, 1), '-.', 'LineWidth', 1.5, 'Color', coltrue);
plot(ground_truth.t, CI_true(:, 2, 2), '-.', 'LineWidth', 1.5, 'Color', coltrue);
h2 = patch([ground_truth.t flip(estimator.results_GTS.t)], [CI_true(:, 2, 1)' flip(CI_true(:, 2, 2))'], ...
           coltrue, 'EdgeAlpha', 0, 'FaceAlpha', .15, 'LineStyle', 'none');

h1 = plot(estimator.results_GTS.t, estimator.results_GTS.population(:, 2), ...
          '-.', 'LineWidth', 2, 'Color', colGTS);
plot(estimator.results_GTS.t, CI_GTS(:, 2, 1), '-.', 'LineWidth', 1.5, 'Color', colGTS);
plot(estimator.results_GTS.t, CI_GTS(:, 2, 2), '-.', 'LineWidth', 1.5, 'Color', colGTS);
h3 = patch([estimator.results_GTS.t flip(estimator.results_GTS.t)], [CI_GTS(:, 2, 1)' flip(CI_GTS(:, 2, 2))'], ...
           colGTS, 'EdgeAlpha', 0, 'FaceAlpha', .15, 'LineStyle', 'none');
       
       
h1 = plot(estimator.results_GMGTS.t_fine, estimator.results_GMGTS.population(:, 2), ...
          ':', 'LineWidth', 2, 'Color', colGMGTS);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 2, 1), ':', 'LineWidth', 1.1, 'Color', colGMGTS);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 2, 2), ':', 'LineWidth', 1.5, 'Color', colGMGTS);
h4 = patch([estimator.results_GMGTS.t_fine flip(estimator.results_GMGTS.t_fine)], [CI_GMGTS(:, 2, 1)' flip(CI_GMGTS(:, 2, 2))'], ...
           colGMGTS, 'EdgeAlpha', 0, 'FaceAlpha', .15, 'LineStyle', 'none');
            
ylim(max(0, ylim))
legend([h(1) h2 h3 h4], {'Data', 'Truth', 'GTS', 'GMGTS'}, 'Location', 'northwest')


%%
filestr = '../../writing/manuscript/figures/maturation.tex';
matlab2tikz(filestr)

fid = fopen(filestr, 'r');
f = fread(fid,'*char')';
f = strrep(f, 'only marks, ', '');
f = strrep(f, 'line width=0.9pt', 'line width=0.5pt, opacity=0.4');
f = strrep(f, 'mark=square, ', 'only marks, mark=square, ');
f = strrep(f, 'mark=triangle, ', 'only marks, mark=triangle, ');
f = strrep(f, 'mark=circle, ', 'only marks, mark=circle, ');
f = strrep(f, 'xlabel style={', 'xlabel style={at={(axis description cs: 0.5, 0.01)}, ');
f = strrep(f, 'ylabel style={', 'ylabel style={at={(axis description cs: 0.03, 0.5)}, ');
f = strrep(f, 'contour/draw color=mycolor1', 'contour/draw color=mycolor1, color=mycolor1');
f = strrep(f, 'contour/draw color=mycolor2', 'contour/draw color=mycolor2, color=mycolor2');
f = strrep(f, 'contour/draw color=mycolor3', 'contour/draw color=mycolor3, color=mycolor3');
f = strrep(f, 'mark size=1.0607pt', 'mark size=2pt');
fclose(fid);
fid = fopen(filestr, 'w');
fwrite(fid, f);
fclose(fid);


