set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

clearvars
close all

load('simulation/bifunctional_measurable.mat'); %10--40-fold

times_GTS(times_GTS > 1000) = mean(times_GTS, 'all');
times_GMGTS_bifunctional = reshape(-times_GMGTS, [], 1);
times_GTS_bifunctional = reshape(times_GTS, [], 1);
xgroupdata_bifunctional = 4*ones(size(times_GMGTS_bifunctional));
cdata_bifunctional = repmat((0:1)', length(times_GMGTS_bifunctional)/2, 1);

load('simulation/maturation_accuracy.mat'); %10--15-fold

times_GMGTS_maturation = reshape(-times_GMGTS, [], 1);
times_GTS_maturation = reshape(times_GTS, [], 1);
xgroupdata_maturation = 3*ones(size(times_GMGTS_maturation));
cdata_maturation = repmat((0:1)', length(times_GMGTS_maturation)/2, 1);

load('simulation/repressilator2.mat'); %6--8-fold

times_GMGTS_repressilator = reshape(-times_GMGTS, [], 1) * .85;
times_GTS_repressilator = reshape(times_GTS, [], 1) * .85;
xgroupdata_repressilator = 2*ones(size(times_GMGTS_repressilator));
cdata_repressilator = repmat((0:1)', length(times_GMGTS_repressilator)/2, 1);


figure('position', [100, 100, 600, 250])
% xgroupdata = repmat(repmat(1:4, 10, 1), 4, 1);
% ydata = -20*ones([size(times_GMGTS, [1 2]) 4]);
% ydata(:, :, 2:end-1) = -times_GMGTS(:, :, :);
% ydata = reshape(permute(ydata, [1 3 2]), [], 4);
% cdata = kron((0:3)', ones(10, 4));
xgroupdata = [xgroupdata_maturation; xgroupdata_bifunctional; xgroupdata_repressilator];
ydata = [times_GMGTS_maturation; times_GMGTS_bifunctional; times_GMGTS_repressilator];
cdata = [cdata_maturation; cdata_bifunctional; cdata_repressilator];

b1 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.', 'BoxWidth', .75);
hold on

% ydata = 40*ones([size(times_GTS, [1 2]) 4]);
% ydata(:, :, 2:end-1) = times_GTS(:, :, :);
% ydata = reshape(permute(ydata, [1 3 2]), [], 4);
% cdata = kron((0:3)', ones(10, 4));
ydata = [times_GTS_maturation; times_GTS_bifunctional; times_GTS_repressilator];

b2 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.', 'BoxWidth', .75);

xline(0, 'k')
xline([-10 10], 'k--')

xlabel("Computing time (sec)")
ylim([0 5])
xlim([-18 70])
xticks(-10:10:70)
yticks([2 3 4])
xticklabels(["10" "0" "10" "20" "30" "40" "50" "60" "70"])
yticklabels(["Repressilator" "Maturation" "Bifunctional TCS"])

b2(1).BoxFaceColor = b1(1).BoxFaceColor;
b2(2).BoxFaceColor = b1(2).BoxFaceColor;
b2(1).MarkerColor = b1(1).MarkerColor;
b2(2).MarkerColor = b1(2).MarkerColor;

% legend([b1(1) b1(2)], ["Full observation" "Partial observation"], 'Location','southeast', 'interpreter', 'latex')
text(-10, .5, "\textbf{GMGTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(10, .5, "\textbf{GTS}", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')



%% SANS SERIF
set(groot,'defaultAxesTickLabelInterpreter','none');  
set(groot,'defaulttextinterpreter','none');
set(groot,'defaultLegendInterpreter','none');


% clearvars
close all

% load('simulation/bifunctional_measurable.mat'); %10--40-fold
% 
% times_GTS(times_GTS > 1000) = mean(times_GTS, 'all');
% times_GMGTS_bifunctional = reshape(-times_GMGTS, [], 1);
% times_GTS_bifunctional = reshape(times_GTS, [], 1);
% xgroupdata_bifunctional = 4*ones(size(times_GMGTS_bifunctional));
% cdata_bifunctional = repmat((0:1)', length(times_GMGTS_bifunctional)/2, 1);
% 
% load('simulation/maturation_accuracy.mat'); %10--15-fold
% 
% times_GMGTS_maturation = reshape(-times_GMGTS, [], 1);
% times_GTS_maturation = reshape(times_GTS, [], 1);
% xgroupdata_maturation = 3*ones(size(times_GMGTS_maturation));
% cdata_maturation = repmat((0:1)', length(times_GMGTS_maturation)/2, 1);
% 
% load('simulation/repressilator2.mat'); %6--8-fold
% 
% times_GMGTS_repressilator = reshape(-times_GMGTS, [], 1) * .85;
% times_GTS_repressilator = reshape(times_GTS, [], 1) * .85;
% xgroupdata_repressilator = 2*ones(size(times_GMGTS_repressilator));
% cdata_repressilator = repmat((0:1)', length(times_GMGTS_repressilator)/2, 1);


figure('position', [100, 100, 400, 250])
% xgroupdata = repmat(repmat(1:4, 10, 1), 4, 1);
% ydata = -20*ones([size(times_GMGTS, [1 2]) 4]);
% ydata(:, :, 2:end-1) = -times_GMGTS(:, :, :);
% ydata = reshape(permute(ydata, [1 3 2]), [], 4);
% cdata = kron((0:3)', ones(10, 4));
xgroupdata = [xgroupdata_maturation; xgroupdata_bifunctional; xgroupdata_repressilator];
ydata = [times_GMGTS_maturation; times_GMGTS_bifunctional; times_GMGTS_repressilator];
cdata = [cdata_maturation; cdata_bifunctional; cdata_repressilator];

b1 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.', 'BoxWidth', .75);
hold on

% ydata = 40*ones([size(times_GTS, [1 2]) 4]);
% ydata(:, :, 2:end-1) = times_GTS(:, :, :);
% ydata = reshape(permute(ydata, [1 3 2]), [], 4);
% cdata = kron((0:3)', ones(10, 4));
ydata = [times_GTS_maturation; times_GTS_bifunctional; times_GTS_repressilator];

b2 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.', 'BoxWidth', .75);

xline(0, 'k')
xline([-10 10], 'k--')

xlabel("Computing time (sec)")
ylim([0 5])
xlim([-18 65])
xticks(-10:10:60)
yticks([2 3 4])
xticklabels(["10" "0" "10" "20" "30" "40" "50" "60"])
yticklabels(["Repressilator" "Maturation" "Bifunctional TCS"])

b2(1).BoxFaceColor = b1(1).BoxFaceColor;
b2(2).BoxFaceColor = b1(2).BoxFaceColor;
b2(1).MarkerColor = b1(1).MarkerColor;
b2(2).MarkerColor = b1(2).MarkerColor;

% legend([b1(1) b1(2)], ["Full observation" "Partial observation"], 'Location','southeast', 'interpreter', 'latex')
text(-10, .5, "GMGTS", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(10, .5, "GTS", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')

%% CONVERT TO LATEX
filestr = '../../writing/manuscript/figures/absolute_times.tex';
% filestr = '../../writing/manuscript/figures/bifunctional_test.tex';
matlab2tikz(filestr)

fid = fopen(filestr, 'r');
f = fread(fid,'*char')';
f = strrep(f, '\begin{axis}[%', '\begin{axis}[scaled ticks=false, tick label style={/pgf/number format/fixed},');
fclose(fid);
fid = fopen(filestr, 'w');
fwrite(fid, f);


