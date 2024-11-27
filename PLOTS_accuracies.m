set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

load('simulation/maturation_accuracy.mat');
% load('simulation/bifunctional_measurable.mat');

close all

first_obs = 2;

figure('Position', [100, 100, 560, 420])
xgroupdata = repmat(repmat(1:6, 10, 1), 6, 1);
ydata = -2*ones([size(hellinger_GMGTS(:, :, :, first_obs), [1 2]) 6]);
ydata(:, :, 2:end-1) = -hellinger_GMGTS(:, :, :, first_obs);
ydata = reshape(permute(ydata, [1 3 2]), [], 6);
cdata = kron((0:5)', ones(10, 6));

b1 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.');
hold on

ydata = 2*ones([size(hellinger_GTS(:, :, :, first_obs), [1 2]) 6]);
ydata(:, :, 2:end-1) = hellinger_GTS(:, :, :, first_obs);
ydata = reshape(permute(ydata, [1 3 2]), [], 6);
cdata = kron((1:6)', ones(10, 6));

b2 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.');

xline(0, 'k')
xline([-.5 .5], 'k--')

ylabel("Noise level")
xlabel("Mismatch (Hellinger distance)")
ylim([.5 6.5])
xlim([-1 1])
xticks(-1:.25:1)
xticklabels(["1" "0.75" "0.5" "0.25" "0" fliplr(["1" "0.75" "0.5" "0.25"])])
yticklabels(["0.1\%", "0.5\%", "1\%", "2\%", "5\%", "10\%"])

b1(5).BoxFaceColor = b1(4).BoxFaceColor;
b1(4).BoxFaceColor = b1(3).BoxFaceColor;
b1(3).BoxFaceColor = b1(2).BoxFaceColor;
b1(2).BoxFaceColor = b1(1).BoxFaceColor;
b1(5).MarkerColor = b1(4).MarkerColor;
b1(4).MarkerColor = b1(3).MarkerColor;
b1(3).MarkerColor = b1(2).MarkerColor;
b1(2).MarkerColor = b1(1).MarkerColor;

legend([b1(5) b1(4) b1(3) b1(2)], ["$T=9$" "$T=13$" "$T=22$" "$T=41$"],'Location','east', 'interpreter', 'latex')
text(-.75, .85, "\textbf{GMGTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(.75, .85, "\textbf{GTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')


%% CONVERT TO LATEX
% % filestr = '../../writing/manuscript/figures/maturation_test.tex';
% filestr = '../../writing/manuscript/figures/bifunctional_test.tex';
% matlab2tikz(filestr)
% 
% fid = fopen(filestr, 'r');
% f = fread(fid,'*char')';
% f = strrep(f, '\begin{axis}[%', '\begin{axis}[scaled ticks=false, tick label style={/pgf/number format/fixed},');
% fclose(fid);
% fid = fopen(filestr, 'w');
% fwrite(fid, f);
% fclose(fid);





%% Wasserstein FULL

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

load('simulation/maturation_accuracy.mat');
% load('simulation/bifunctional_measurable.mat');

close all

first_obs = 1;

figure('Position', [100, 100, 350, 250])
xgroupdata = repmat(repmat(1:6, 10, 1), 6, 1);
ydata = -2*ones([size(wasserstein_GMGTS(:, :, :, first_obs), [1 2]) 6]);
ydata(:, :, 2:end-1) = -wasserstein_GMGTS(:, :, :, first_obs);
ydata = reshape(permute(ydata, [1 3 2]), [], 6);
cdata = kron((0:5)', ones(10, 6));

b1 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.');
hold on

ydata = 2*ones([size(wasserstein_GTS(:, :, :, first_obs), [1 2]) 6]);
ydata(:, :, 2:end-1) = wasserstein_GTS(:, :, :, first_obs);
ydata = reshape(permute(ydata, [1 3 2]), [], 6);
cdata = kron((1:6)', ones(10, 6));

b2 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.');

xline(0, 'k')

ylabel("Noise level")
xlabel("Mismatch (Wasserstein distance)")
ylim([.5 6.5])

% xlim([-.32 .15])
% xticks(-.3:.1:.1)
% xticklabels(["0.3" "0.2" "0.1" "0" "0.1"])
% xline([-.05 .05], 'k--')
% text(-.125, .85, "\textbf{GMGTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
% text(.1, .85, "\textbf{GTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')

xlim([-.006 .006])
xticks(-.006:.002:.006)
xticklabels(["0.006" "0.004" "0.002" "0" fliplr(["0.006" "0.004" "0.002"])])
xline([-.002 .002], 'k--')
text(-.004, .85, "\textbf{GMGTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(.004, .85, "\textbf{GTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')

yticklabels(["0.1\%", "0.5\%", "1\%", "2\%", "5\%", "10\%"])

b1(5).BoxFaceColor = b1(4).BoxFaceColor;
b1(4).BoxFaceColor = b1(3).BoxFaceColor;
b1(3).BoxFaceColor = b1(2).BoxFaceColor;
b1(2).BoxFaceColor = b1(1).BoxFaceColor;
b1(5).MarkerColor = b1(4).MarkerColor;
b1(4).MarkerColor = b1(3).MarkerColor;
b1(3).MarkerColor = b1(2).MarkerColor;
b1(2).MarkerColor = b1(1).MarkerColor;

% legend([b1(5) b1(4) b1(3) b1(2)], ["$T=9$" "$T=13$" "$T=22$" "$T=41$"], 'Location', 'west', 'interpreter', 'latex')



%% Wasserstein PARTIAL

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

load('simulation/maturation_accuracy.mat');
% load('simulation/bifunctional_measurable.mat');

close all

first_obs = 2;

figure('Position', [100, 100, 350, 250])
xgroupdata = repmat(repmat(1:6, 10, 1), 6, 1);
ydata = -2*ones([size(wasserstein_GMGTS(:, :, :, first_obs), [1 2]) 6]);
ydata(:, :, 2:end-1) = -wasserstein_GMGTS(:, :, :, first_obs);
ydata = reshape(permute(ydata, [1 3 2]), [], 6);
cdata = kron((0:5)', ones(10, 6));

b1 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.');
hold on

ydata = 2*ones([size(wasserstein_GTS(:, :, :, first_obs), [1 2]) 6]);
ydata(:, :, 2:end-1) = wasserstein_GTS(:, :, :, first_obs);
ydata = reshape(permute(ydata, [1 3 2]), [], 6);
cdata = kron((1:6)', ones(10, 6));

b2 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.');

xline(0, 'k')



ylabel("Noise level")
xlabel("Mismatch (Wasserstein distance)")
ylim([.5 6.5])

% xlim([-.35 .35])
% xticks(-.3:.1:.3)
% xticklabels(["0.3" "0.2" "0.1" "0" fliplr(["0.3" "0.2" "0.1"])])
% xline([-.1 .1], 'k--')
% text(-.25, .85, "\textbf{GMGTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
% text(.25, .85, "\textbf{GTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')

xlim([-.015 .01])
xticks(-.012:.004:.01)
xticklabels(["0.0012" "0.008" "0.004" "0" fliplr(["0.008" "0.004"])])
xline([-.004 .004], 'k--')
text(-.009, .85, "\textbf{GMGTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(.008, .85, "\textbf{GTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')


yticklabels(["0.1\%", "0.5\%", "1\%", "2\%", "5\%", "10\%"])

b1(5).BoxFaceColor = b1(4).BoxFaceColor;
b1(4).BoxFaceColor = b1(3).BoxFaceColor;
b1(3).BoxFaceColor = b1(2).BoxFaceColor;
b1(2).BoxFaceColor = b1(1).BoxFaceColor;
b1(5).MarkerColor = b1(4).MarkerColor;
b1(4).MarkerColor = b1(3).MarkerColor;
b1(3).MarkerColor = b1(2).MarkerColor;
b1(2).MarkerColor = b1(1).MarkerColor;

% legend([b1(5) b1(4) b1(3) b1(2)], ["$T=9$" "$T=13$" "$T=22$" "$T=41$"], 'Location', 'east', 'interpreter', 'latex')



%% Wasserstein PARTIAL

set(groot,'defaultAxesTickLabelInterpreter','none');  
set(groot,'defaulttextinterpreter','none');
set(groot,'defaultLegendInterpreter','none');

load('simulation/maturation_accuracy.mat');
% load('simulation/bifunctional_measurable.mat');

close all

first_obs = 2;

% figure('Position', [100, 100, 350, 250])
figure('Position', [100, 100, 450, 350])
xgroupdata = repmat(repmat(1:6, 10, 1), 6, 1);
ydata = -2*ones([size(wasserstein_GMGTS(:, :, :, first_obs), [1 2]) 6]);
ydata(:, :, 2:end-1) = -wasserstein_GMGTS(:, :, :, first_obs);
ydata = reshape(permute(ydata, [1 3 2]), [], 6);
cdata = kron((0:5)', ones(10, 6));

b1 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.');
hold on

ydata = 2*ones([size(wasserstein_GTS(:, :, :, first_obs), [1 2]) 6]);
ydata(:, :, 2:end-1) = wasserstein_GTS(:, :, :, first_obs);
ydata = reshape(permute(ydata, [1 3 2]), [], 6);
cdata = kron((1:6)', ones(10, 6));

b2 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.');

xline(0, 'k')



ylabel("Noise level")
xlabel("Mismatch (Wasserstein distance)")
ylim([.5 6.5])

% xlim([-.35 .35])
% xticks(-.3:.1:.3)
% xticklabels(["0.3" "0.2" "0.1" "0" fliplr(["0.3" "0.2" "0.1"])])
% xline([-.1 .1], 'k--')
% text(-.25, .85, "\textbf{GMGTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
% text(.25, .85, "\textbf{GTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')

xlim([-.015 .01])
xticks(-.012:.004:.01)
xticklabels(["0.0012" "0.008" "0.004" "0" fliplr(["0.008" "0.004"])])
xline([-.004 .004], 'k--')
text(-.009, .85, "GMGTS", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(.008, .85, "GTS", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')


yticklabels(["0.1%", "0.5%", "1%", "2%", "5%", "10%"])

b1(5).BoxFaceColor = b1(4).BoxFaceColor;
b1(4).BoxFaceColor = b1(3).BoxFaceColor;
b1(3).BoxFaceColor = b1(2).BoxFaceColor;
b1(2).BoxFaceColor = b1(1).BoxFaceColor;
b1(5).MarkerColor = b1(4).MarkerColor;
b1(4).MarkerColor = b1(3).MarkerColor;
b1(3).MarkerColor = b1(2).MarkerColor;
b1(2).MarkerColor = b1(1).MarkerColor;

legend([b1(5) b1(4) b1(3) b1(2)], ["$T=9$" "$T=13$" "$T=22$" "$T=41$"], 'Location', 'west', 'interpreter', 'latex')