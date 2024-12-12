%% Wasserstein FULL STACKED

set(groot,'defaultAxesTickLabelInterpreter','none');  
set(groot,'defaulttextinterpreter','none');
set(groot,'defaultLegendInterpreter','none');

load('simulation/maturation_accuracy2.mat');
% load('simulation/bifunctional_accuracy2.mat');

close all

first_obs = 1;

figure('Position', [100, 100, 350, 500])
xgroupdata = repmat(repmat(1:12, 10, 1), 6, 1);
ydata = -2*ones([size(ws_GMGTS(:, :, :, first_obs), [1 2]) 6]);
ydata(:, :, 2:end-1) = ws_GMGTS(:, :, :, first_obs);
ydata = reshape(permute(ydata, [1 3 2]), [], 6);
cdata = kron((0:5)', ones(10, 6));

ydata2 = 2*ones([size(ws_GTS(:, :, :, first_obs), [1 2]) 6]);
ydata2(:, :, 2:end-1) = ws_GTS(:, :, :, first_obs);
ydata2 = reshape(permute(ydata2, [1 3 2]), [], 6);
cdata2 = kron((0:5)', ones(10, 6));

ydata = [ydata ydata2];
cdata = [cdata cdata2];

b1 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.');

xline(0, 'k')
yline(6.5, 'k')

ylabel("Noise level")
xlabel("Mismatch (Wasserstein distance)")
ylim([.5 12.5])

xlim([0 .6])
xline(.1:.1:.3, '-', Color=[.5 .5 .5])
text(.5, .85, "GMGTS", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(.5, 6.85, "GTS", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
yticks(1:12)
yticklabels(repmat(["0.1%" "0.5%" "1%" "2%" "5%" "10%"], 1, 2))

b1(5).BoxFaceColor = b1(4).BoxFaceColor;
b1(4).BoxFaceColor = b1(3).BoxFaceColor;
b1(3).BoxFaceColor = b1(2).BoxFaceColor;
b1(2).BoxFaceColor = b1(1).BoxFaceColor;
b1(5).MarkerColor = b1(4).MarkerColor;
b1(4).MarkerColor = b1(3).MarkerColor;
b1(3).MarkerColor = b1(2).MarkerColor;
b1(2).MarkerColor = b1(1).MarkerColor;

legend([b1(5) b1(4) b1(3) b1(2)], ["$\mathsf{T=9}$" "$\mathsf{T=13}$" "$\mathsf{T=22}$" "$\mathsf{T=41}$"], 'Location', 'northeast', 'interpreter', 'latex')



%% Wasserstein PARTIAL STACKED

set(groot,'defaultAxesTickLabelInterpreter','none');  
set(groot,'defaulttextinterpreter','none');
set(groot,'defaultLegendInterpreter','none');

load('simulation/maturation_accuracy.mat');
% load('simulation/bifunctional_accuracy2.mat');

% close all

first_obs = 2;

figure('Position', [100, 100, 350, 500])
xgroupdata = repmat(repmat(1:12, 10, 1), 6, 1);
ydata = -2*ones([size(ws_GMGTS(:, :, :, first_obs), [1 2]) 6]);
ydata(:, :, 2:end-1) = ws_GMGTS(:, :, :, first_obs);
ydata = reshape(permute(ydata, [1 3 2]), [], 6);
cdata = kron((0:5)', ones(10, 6));

ydata2 = 2*ones([size(ws_GTS(:, :, :, first_obs), [1 2]) 6]);
ydata2(:, :, 2:end-1) = ws_GTS(:, :, :, first_obs);
ydata2 = reshape(permute(ydata2, [1 3 2]), [], 6);
cdata2 = kron((0:5)', ones(10, 6));

ydata = [ydata ydata2];
cdata = [cdata cdata2];

b1 = boxchart(xgroupdata(:), ydata(:), 'GroupByColor', cdata(:), 'Orientation', 'horizontal', 'MarkerStyle', '.');

xline(0, 'k')
yline(6.5, 'k')

ylabel("Noise level")
xlabel("Mismatch (Wasserstein distance)")
ylim([.5 12.5])

xlim([0 .8])
xline(.1:.1:.3, '-', Color=[.5 .5 .5])
text(.8*5/6, .85, "GMGTS", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(.8*5/6, 6.85, "GTS", 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
yticks(1:12)
yticklabels(repmat(["0.1%" "0.5%" "1%" "2%" "5%" "10%"], 1, 2))

b1(5).BoxFaceColor = b1(4).BoxFaceColor;
b1(4).BoxFaceColor = b1(3).BoxFaceColor;
b1(3).BoxFaceColor = b1(2).BoxFaceColor;
b1(2).BoxFaceColor = b1(1).BoxFaceColor;
b1(5).MarkerColor = b1(4).MarkerColor;
b1(4).MarkerColor = b1(3).MarkerColor;
b1(3).MarkerColor = b1(2).MarkerColor;
b1(2).MarkerColor = b1(1).MarkerColor;

% legend([b1(5) b1(4) b1(3) b1(2)], ["$T=9$" "$T=13$" "$T=22$" "$T=41$"], 'Location', 'northeast', 'interpreter', 'latex')