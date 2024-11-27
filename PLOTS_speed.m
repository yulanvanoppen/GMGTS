set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

load('simulation/genlotka_speed.mat');
% load('simulation/bifunctional_measurable2.mat');

close all

figure('position', [100, 100, 600, 250])
% xgroupdata = repmat(repmat(1:4, 10, 1), 4, 1);
% ydata = -20*ones([size(times_GMGTS, [1 2]) 4]);
% ydata(:, :, 2:end-1) = -times_GMGTS(:, :, :);
% ydata = reshape(permute(ydata, [1 3 2]), [], 4);
% cdata = kron((0:3)', ones(10, 4));
xgroupdata = repmat(repmat([2 3 5 9], 10, 1), 2, 1);
ydata = -times_GMGTS(:, :, :);
ydata = reshape(permute(ydata, [1 3 2]), [], 4);
cdata = kron((0:1)', ones(10, 4));

b1 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.', 'BoxWidth', 1);
hold on
l1 = plot(-permute(median(times_GMGTS), [2 3 1]), [[2 3 5 9]-.25; [2 3 5 9]+.25]', ':', 'LineWidth', 1.5);
set(l1, {'color'}, num2cell([b1(1).MarkerColor; b1(2).MarkerColor], 2))

% ydata = 40*ones([size(times_GTS, [1 2]) 4]);
% ydata(:, :, 2:end-1) = times_GTS(:, :, :);
% ydata = reshape(permute(ydata, [1 3 2]), [], 4);
% cdata = kron((0:3)', ones(10, 4));
ydata = times_GTS(:, :, :);
ydata = reshape(permute(ydata, [1 3 2]), [], 2);
cdata = kron((0:1)', ones(10, 4));

b2 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.', 'BoxWidth', 1);
l2 = plot(permute(median(times_GTS), [2 3 1]), [[2 3 5 9]-.25; [2 3 5 9]+.25]', ':', 'LineWidth', 1.5);
set(l2, {'color'}, num2cell([b1(1).MarkerColor; b1(2).MarkerColor], 2))

xline(0, 'k')
xline([-10 10], 'k--')

xlabel("Computing time (sec)")
ylim([0 10])
xlim([-15 45])
xticks(-15:5:45)
yticks([2 3 5 9])
xticklabels(["15" "10" "5" "0" "5" "10" "15" "20" "25" "30" "35" "40" "45"])
yticklabels(["$P=1$" "$P=2$" "$P=4$" "$P=8$"])

b2(1).BoxFaceColor = b1(1).BoxFaceColor;
b2(2).BoxFaceColor = b1(2).BoxFaceColor;
b2(1).MarkerColor = b1(1).MarkerColor;
b2(2).MarkerColor = b1(2).MarkerColor;

legend([b1(2) b1(1)], ["Partial observation" "Full observation"], 'Location','southeast', 'interpreter', 'latex')
text(-5, .5, "\textbf{GMGTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(5, .5, "\textbf{GTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')


%% CONVERT TO LATEX
% filestr = '../../writing/manuscript/figures/maturation_test.tex';
% % filestr = '../../writing/manuscript/figures/bifunctional_test.tex';
% matlab2tikz(filestr)
% 
% fid = fopen(filestr, 'r');
% f = fread(fid,'*char')';
% f = strrep(f, '\begin{axis}[%', '\begin{axis}[scaled ticks=false, tick label style={/pgf/number format/fixed},');
% fclose(fid);
% fid = fopen(filestr, 'w');
% fwrite(fid, f);



set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

load('simulation/genlotka_speed.mat');
% load('simulation/bifunctional_measurable2.mat');

close all

figure('position', [100, 100, 600, 250])
% xgroupdata = repmat(repmat(1:4, 10, 1), 4, 1);
% ydata = -20*ones([size(times_GMGTS, [1 2]) 4]);
% ydata(:, :, 2:end-1) = -times_GMGTS(:, :, :);
% ydata = reshape(permute(ydata, [1 3 2]), [], 4);
% cdata = kron((0:3)', ones(10, 4));
xgroupdata = repmat(repmat([2 3 5 9], 10, 1), 2, 1);
ydata = -accuracies_GMGTS(:, :, :);
ydata = reshape(permute(ydata, [1 3 2]), [], 2);
cdata = kron((0:1)', ones(10, 4));

b1 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.', 'BoxWidth', 1);
hold on
l1 = plot(-permute(median(accuracies_GMGTS), [2 3 1]), [[2 3 5 9]-.25; [2 3 5 9]+.25]', ':', 'LineWidth', 1.5);
set(l1, {'color'}, num2cell([b1(1).MarkerColor; b1(2).MarkerColor], 2))

% ydata = 40*ones([size(times_GTS, [1 2]) 4]);
% ydata(:, :, 2:end-1) = times_GTS(:, :, :);
% ydata = reshape(permute(ydata, [1 3 2]), [], 4);
% cdata = kron((0:3)', ones(10, 4));
ydata = accuracies_GTS(:, :, :);
ydata = reshape(permute(ydata, [1 3 2]), [], 2);
cdata = kron((0:1)', ones(10, 4));

b2 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.', 'BoxWidth', 1);
l2 = plot(permute(median(accuracies_GTS), [2 3 1]), [[2 3 5 9]-.25; [2 3 5 9]+.25]', ':', 'LineWidth', 1.5);
set(l2, {'color'}, num2cell([b1(1).MarkerColor; b1(2).MarkerColor], 2))

xline(0, 'k')
xline([-.0002 .0002], 'k--')

xlabel("Mismatch (Wasserstein distance)")
ylim([0 10])
xlim([-0.001 0.0005])
xticks(-0.001:0.0002:0.0005)
yticks([2 3 5 9])
xticklabels(["1E-3" "8E-4" "6E-4" "4E-4" "2E-4" "0" "2E-4" "4E-4"])
yticklabels(["$P=1$" "$P=2$" "$P=4$" "$P=8$"])

b2(1).BoxFaceColor = b1(1).BoxFaceColor;
b2(2).BoxFaceColor = b1(2).BoxFaceColor;
b2(1).MarkerColor = b1(1).MarkerColor;
b2(2).MarkerColor = b1(2).MarkerColor;

% legend([b1(1) b1(2)], ["Full observation" "Partial observation"], 'Location','east', 'interpreter', 'latex')
text(-.0004, .5, "\textbf{GMGTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(.0004, .5, "\textbf{GTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')



%% SANS SERIF
set(groot,'defaultAxesTickLabelInterpreter','tex');  
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

load('simulation/genlotka_speed.mat');
% load('simulation/bifunctional_measurable2.mat');

close all

figure('position', [100, 100, 400, 250])
% xgroupdata = repmat(repmat(1:4, 10, 1), 4, 1);
% ydata = -20*ones([size(times_GMGTS, [1 2]) 4]);
% ydata(:, :, 2:end-1) = -times_GMGTS(:, :, :);
% ydata = reshape(permute(ydata, [1 3 2]), [], 4);
% cdata = kron((0:3)', ones(10, 4));
xgroupdata = repmat(repmat([2 3 5 9], 10, 1), 2, 1);
ydata = -times_GMGTS(:, :, :);
ydata = reshape(permute(ydata, [1 3 2]), [], 4);
cdata = kron((0:1)', ones(10, 4));

b1 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.', 'BoxWidth', 1);
hold on
l1 = plot(-permute(median(times_GMGTS), [2 3 1]), [[2 3 5 9]-.25; [2 3 5 9]+.25]', ':', 'LineWidth', 1.5);
set(l1, {'color'}, num2cell([b1(1).MarkerColor; b1(2).MarkerColor], 2))

% ydata = 40*ones([size(times_GTS, [1 2]) 4]);
% ydata(:, :, 2:end-1) = times_GTS(:, :, :);
% ydata = reshape(permute(ydata, [1 3 2]), [], 4);
% cdata = kron((0:3)', ones(10, 4));
ydata = times_GTS(:, :, :);
ydata = reshape(permute(ydata, [1 3 2]), [], 2);
cdata = kron((0:1)', ones(10, 4));

b2 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.', 'BoxWidth', 1);
l2 = plot(permute(median(times_GTS), [2 3 1]), [[2 3 5 9]-.25; [2 3 5 9]+.25]', ':', 'LineWidth', 1.5);
set(l2, {'color'}, num2cell([b1(1).MarkerColor; b1(2).MarkerColor], 2))

xline(0, 'k')
xline([-10 10], 'k--')

xlabel("Computing time (sec)")
ylim([0 10])
xlim([-15 45])
xticks(-15:5:45)
yticks([2 3 5 9])
xticklabels(["15" "10" "5" "0" "5" "10" "15" "20" "25" "30" "35" "40" "45"])
yticklabels(["{\it P}=1" "{\it P}=2" "{\it P}=4" "{\it P}=8"])

b2(1).BoxFaceColor = b1(1).BoxFaceColor;
b2(2).BoxFaceColor = b1(2).BoxFaceColor;
b2(1).MarkerColor = b1(1).MarkerColor;
b2(2).MarkerColor = b1(2).MarkerColor;

legend([b1(2) b1(1)], ["Partial obs." "Full obs."], 'Location','southeast')
text(-7.5, .75, "GMGTS", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(7.5, .75, "GTS", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')


%% CONVERT TO LATEX
% filestr = '../../writing/manuscript/figures/maturation_test.tex';
% % filestr = '../../writing/manuscript/figures/bifunctional_test.tex';
% matlab2tikz(filestr)
% 
% fid = fopen(filestr, 'r');
% f = fread(fid,'*char')';
% f = strrep(f, '\begin{axis}[%', '\begin{axis}[scaled ticks=false, tick label style={/pgf/number format/fixed},');
% fclose(fid);
% fid = fopen(filestr, 'w');
% fwrite(fid, f);
% fclose(fid);