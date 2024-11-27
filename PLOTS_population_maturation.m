clearvars

% files = dir(fullfile('simulation', '*.mat'));
% files = {files.name};

load(['simulation/maturation'])


%%
close all

figure('position', [70, 10, 800, 600])
tiledlayout(2, 3);

legends = [];
parameter_names = system.parameters_variable;

labels = {'Data', 'GMGTS', 'GTS', 'Truth'};

nexttile(1)
title(system.states{1})
xlabel('Time (min)')
ylabel('Measurement (A.U.)')
hold on

CI_GMGTS = quantile(estimator.results_GMGTS.fitted2, [.05 .95], 3);
%             switch marker
%                 case ':', col = [0, 0.4470, 0.7410];
%                 case '--', col = [0.8500, 0.3250, 0.0980];
%                 case '-.', col = [0.4940, 0.1840, 0.5560];
%             end
            
h1 = plot(estimator.results_GMGTS.t_fine, estimator.results_GMGTS.population(:, 1), ...
          ':', 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 1, 1), ':', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 1, 2), ':', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
h2 = patch([estimator.results_GMGTS.t_fine flip(estimator.results_GMGTS.t_fine)], [CI_GMGTS(:, 1, 1)' flip(CI_GMGTS(:, 1, 2))'], ...
           [0, 0.4470, 0.7410], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');
       

CI_GTS = quantile(estimator.results_GTS.fitted2, [.05 .95], 3);
%             switch marker
%                 case ':', col = [0, 0.4470, 0.7410];
%                 case '--', col = [0.8500, 0.3250, 0.0980];
%                 case '-.', col = [0.4940, 0.1840, 0.5560];
%             end
            
h1 = plot(estimator.results_GTS.t, estimator.results_GTS.population(:, 1), ...
          '-.', 'LineWidth', 2, 'Color', [0.4940, 0.1840, 0.5560]);
plot(estimator.results_GTS.t, CI_GTS(:, 1, 1), '-.', 'LineWidth', 1.5, 'Color', [0.4940, 0.1840, 0.5560]);
plot(estimator.results_GTS.t, CI_GTS(:, 1, 2), '-.', 'LineWidth', 1.5, 'Color', [0.4940, 0.1840, 0.5560]);
h2 = patch([estimator.results_GTS.t flip(estimator.results_GTS.t)], [CI_GTS(:, 1, 1)' flip(CI_GTS(:, 1, 2))'], ...
           [0.4940, 0.1840, 0.5560], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');
       

CI_true = quantile(ground_truth.original, [.05 .95], 3);
%             switch marker
%                 case ':', col = [0, 0.4470, 0.7410];
%                 case '--', col = [0.8500, 0.3250, 0.0980];
%                 case '-.', col = [0.4940, 0.1840, 0.5560];
%             end
            
h1 = plot(data.t, ground_truth.moriginal(:, 1), ...
          ':', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
plot(data.t, CI_true(:, 1, 1), '--', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
plot(data.t, CI_true(:, 1, 2), '--', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
h2 = patch([data.t flip(data.t)], [CI_true(:, 1, 1)' flip(CI_true(:, 1, 2))'], ...
           [0.8500 0.3250 0.0980], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');



nexttile(4)
title(system.states{5})
xlabel('Time (min)')
ylabel('Measurement (A.U.)')
hold on
            
h1 = plot(estimator.results_GMGTS.t_fine, estimator.results_GMGTS.population(:, 5), ...
          ':', 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 5, 1), ':', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 5, 2), ':', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
h2 = patch([estimator.results_GMGTS.t_fine flip(estimator.results_GMGTS.t_fine)], [CI_GMGTS(:, 5, 1)' flip(CI_GMGTS(:, 5, 2))'], ...
           [0, 0.4470, 0.7410], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');
       

h1 = plot(estimator.results_GTS.t, estimator.results_GTS.population(:, 5), ...
          '-.', 'LineWidth', 2, 'Color', [0.4940, 0.1840, 0.5560]);
plot(estimator.results_GTS.t, CI_GTS(:, 5, 1), '-.', 'LineWidth', 1.5, 'Color', [0.4940, 0.1840, 0.5560]);
plot(estimator.results_GTS.t, CI_GTS(:, 5, 2), '-.', 'LineWidth', 1.5, 'Color', [0.4940, 0.1840, 0.5560]);
h2 = patch([estimator.results_GTS.t flip(estimator.results_GTS.t)], [CI_GTS(:, 5, 1)' flip(CI_GTS(:, 5, 2))'], ...
           [0.4940, 0.1840, 0.5560], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');
       

h1 = plot(data.t, ground_truth.moriginal(:, 5), ...
          ':', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
plot(data.t, CI_true(:, 5, 1), '--', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
plot(data.t, CI_true(:, 5, 2), '--', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
h2 = patch([data.t flip(data.t)], [CI_true(:, 5, 1)' flip(CI_true(:, 5, 2))'], ...
           [0.8500 0.3250 0.0980], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');




nexttile([2 2])
title(system.states{6})
xlabel('Time (min)')
ylabel('Measurement (A.U.)')
hold on

h = plot(data.t, reshape(data.traces(:, 1, 1:data.N), [], data.N), '-');
set(h, {'color'}, num2cell(parula(data.N), 2));

h0 = plot(estimator.results_GMGTS.t_fine, estimator.results_GMGTS.population(:, 6), ...
          ':', 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 6, 1), ':', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
plot(estimator.results_GMGTS.t_fine, CI_GMGTS(:, 6, 2), ':', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
h1 = patch([estimator.results_GMGTS.t_fine flip(estimator.results_GMGTS.t_fine)], [CI_GMGTS(:, 6, 1)' flip(CI_GMGTS(:, 6, 2))'], ...
           [0, 0.4470, 0.7410], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');
       

h0 = plot(estimator.results_GTS.t, estimator.results_GTS.population(:, 6), ...
          '-.', 'LineWidth', 2, 'Color', [0.4940, 0.1840, 0.5560]);
plot(estimator.results_GTS.t, CI_GTS(:, 6, 1), '-.', 'LineWidth', 1.5, 'Color', [0.4940, 0.1840, 0.5560]);
plot(estimator.results_GTS.t, CI_GTS(:, 6, 2), '-.', 'LineWidth', 1.5, 'Color', [0.4940, 0.1840, 0.5560]);
h2 = patch([estimator.results_GTS.t flip(estimator.results_GTS.t)], [CI_GTS(:, 6, 1)' flip(CI_GTS(:, 6, 2))'], ...
           [0.4940, 0.1840, 0.5560], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');
       

h0 = plot(data.t, ground_truth.moriginal(:, 6), ...
          ':', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
plot(data.t, CI_true(:, 6, 1), '--', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
plot(data.t, CI_true(:, 6, 2), '--', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
h3 = patch([data.t flip(data.t)], [CI_true(:, 6, 1)' flip(CI_true(:, 6, 2))'], ...
           [0.8500 0.3250 0.0980], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');
       
legend([h(1) h1 h2 h3], labels, 'Location', 'northwest')




% h = plot(data.t, reshape(data.traces(:, 1, 1:data.N), [], data.N), '-');
% set(h, {'color'}, num2cell(parula(data.N), 2));
%             
% h0 = plot(estimator.results_GMGTS.t, estimator.results_GMGTS.population(:, 6), ...
%           ':', 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
% plot(estimator.results_GMGTS.t, CI_GMGTS(:, 6, 1), ':', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
% plot(estimator.results_GMGTS.t, CI_GMGTS(:, 6, 2), ':', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
% h1 = patch([estimator.results_GMGTS.t flip(estimator.results_GMGTS.t)], [CI_GMGTS(:, 6, 1)' flip(CI_GMGTS(:, 6, 2))'], ...
%            [0, 0.4470, 0.7410], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');
%             
% h0 = plot(estimator.results_GTS.t, estimator.results_GTS.population(:, 6), ...
%           '-.', 'LineWidth', 2, 'Color', [0.4940, 0.1840, 0.5560]);
% plot(estimator.results_GTS.t, CI_GTS(:, 6, 1), '-.', 'LineWidth', 1.5, 'Color', [0.4940, 0.1840, 0.5560]);
% plot(estimator.results_GTS.t, CI_GTS(:, 6, 2), '-.', 'LineWidth', 1.5, 'Color', [0.4940, 0.1840, 0.5560]);
% h2 = patch([estimator.results_GTS.t flip(estimator.results_GTS.t)], [CI_GTS(:, 6, 1)' flip(CI_GTS(:, 6, 2))'], ...
%            [0.4940, 0.1840, 0.5560], 'EdgeAlpha', 0, 'FaceAlpha', .3, 'LineStyle', 'none');
%        
% legend([h(1) h1 h2], labels, 'Location', 'northwest')


%%
filestr = '../../writing/manuscript/figures/maturation_simulation_population.tex';
matlab2tikz(filestr)

fid = fopen(filestr, 'r');
f = fread(fid,'*char')';
f = strrep(f, 'only marks, ', '');
fclose(fid);
fid = fopen(filestr, 'w');
fwrite(fid, f);
fclose(fid);




