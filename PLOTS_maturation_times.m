clearvars

folder = 'experimental';
files = dir(fullfile(folder, '*.mat'));
files = {files.name};
twostep = logical([0 1 0 1 0 1 1 0 0 0 1 1]);

nfiles = length(files);

lmeans = zeros(nfiles, 1);
lsds = zeros(nfiles, 1);

for index = 1:nfiles
    load([folder '/' files{index}])
    lmeans(index) = estimator.results_GMGTS.b_est(end);
    lsds(index) = sqrt(estimator.results_GMGTS.D_est(end));
end

colors_hex = ["#00ccff"; % CFP
              "#f70000"; % mCherry
              "#d6ff00"; % mCitrine
              "#d70000"; % mKate2
              "#19ff00"; % mNeonGreen
              "#ff2f00"; % mScarletI
              "#00adff"; % mTurquoise2
              "#cfff00"; % mVenus
              "#12ff00"; % pH-tdGFP
              "#00ff00"; % sfGFP
              "#00ffe5"; % mTFP1
              "#ff5e00"]; % tdTomato]

colors_hex2 = arrayfun(@(x) eraseBetween(x, 1, 1), colors_hex);
colors_triplets = arrayfun(@hex2rgb, colors_hex2, 'UniformOutput', false);

means = exp(lmeans + lsds.^2/2);
factor = 1.67835 * twostep' + log(2) * ~twostep';
[ht_means, order] = sort(means ./ factor, 'ascend');
means = means(order);
lmeans = lmeans(order);
lsds = lsds(order);

twostep = twostep(order);
colors_triplets = colors_triplets(order);

gray = 0.3;



%%
close all

% km_range = [means(1)-.5*sds(1), means(end)+10*sds(end)];
t_points = 601;
t_range = [1e-1 300];
t_grid = linspace(t_range(1), t_range(2), t_points);
mtime_grid = log(2) * t_grid;
mtime2_grid = 1.67835 * t_grid;

densities = zeros(t_points, nfiles);

for index = 1:nfiles
    densities(:, index) = lognpdf(t_grid, -lmeans(index), lsds(index));
end

max_densities =  max(densities);
densities_normalized = densities ./ max_densities;
densities_corrected = densities_normalized * .9;

figure('position', [100, 100, 500, 300])
hold on
for index = nfiles:-1:1
    col = colors_triplets{index};
    col = (1-gray) * col + gray * [1 1 1];
    
    if twostep(index), tgrid = mtime2_grid; else, tgrid = mtime_grid; end
    plot(tgrid, densities_corrected(:, index)+index/2, 'Color', col*.8, 'LineWidth', 1.25)
    patch([tgrid, flip(tgrid)], [densities_corrected(:, index)'+index/2, index/2*ones(1, t_points)], ...
          col*.9, 'EdgeAlpha', 0, 'FaceAlpha', .5)
      
    factor = 1.67835 * twostep(index) + log(2) * ~twostep(index);
    plot(factor/means(index), ...
        .9*lognpdf(1/means(index), -lmeans(index), lsds(index)) / max_densities(index) + index/2, ...
         '.', 'MarkerSize', 15, 'Color', col*.8);
end
hold off

names2 = string(files(order));
for index = 1:nfiles
    names2(index) = regexprep(names2(index), "EL222_AQTrip[-_]", "");
    names2(index) = regexprep(names2(index), "_.+", "");
end

grid on
xticks([.5 1 2 5 10 20 50 100 200 300])
yticks([(1:nfiles)/2 nfiles/2+1])
yticklabels([names2 ""])
xlabel('Maturation time (min)')
xlim([1 200])
ylim([0 nfiles/2 + 1])
set(gca, 'XScale', 'log')



%%
filestr = 'figures/maturation_times_sf.tex';
matlab2tikz(filestr)

fid = fopen(filestr, 'r');
f = fread(fid,'*char')';
f = strrep(f, 'xminorticks=true,', ...
           'xminorticks=true, minor xtick={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200},');
f = strrep(f, 'xmode=log,', 'xmode=log, xticklabels={{.5}, {1}, {2}, {5}, {10}, {20}, {50}, {100}, {200}},');
fclose(fid);
fid = fopen(filestr, 'w');
fwrite(fid, f);
fclose(fid);



%%
times_paolo = [.084 .034 .058 .015 .058 .058 .036 .05 .056 .1 .021 .016];
twostep_indices = [2 4 6 7 11 12];
times_paolo = log(2)./times_paolo;
times_paolo(twostep_indices) = times_paolo(twostep_indices) * 1.67835 / log(2);
times_paolo = round(times_paolo(order), 1);

modes_quartiles_times = zeros(nfiles, 5);
for index = 1:nfiles
    factor = 1.67835 * twostep(nfiles-index+1) + log(2) * ~twostep(nfiles-index+1);
    modes_quartiles_times(index, 1) = factor/means(nfiles-index+1);
    modes_quartiles_times(index, 2:3) = factor./exp(icdf('Normal', [.75 .25], lmeans(nfiles-index+1), lsds(nfiles-index+1)));
    modes_quartiles_times(index, 4) = diff(modes_quartiles_times(index, 2:3));
    modes_quartiles_times(index, 5) = .75 * modes_quartiles_times(index, 4) / modes_quartiles_times(index, 1);
end
modes_quartiles_times(:, 1:4) = round(modes_quartiles_times(:, 1:4), 1);
modes_quartiles_times(:, 5) = round(modes_quartiles_times(:, 5), 3)
dat = table(flip(names2'), flip(times_paolo'), modes_quartiles_times(:, [1:3 5]), ...
            'VariableNames', {'Protein', 'Paolo', 'GMGTS_q25_q75_IQR_RQR'})
        
protein = flip(names2');
paolos = flip(times_paolo');
modes = modes_quartiles_times(:, 1);
q25 = modes_quartiles_times(:, 2);
q75 = modes_quartiles_times(:, 3);
times = modes_quartiles_times(:, 5);
table2latex(protein, paolos, modes, q25, q75, times)
