clearvars
close all
rng('default')

files = dir(fullfile('estimates_final3', '*.mat'));
files = {files.name}';
nfiles = length(files);

nrep = 10;
loglik_fixed = zeros(nfiles, nrep);
loglik_final = zeros(nfiles, nrep);
loglik_final_GTS = zeros(nfiles, nrep);

for index = 10
    fprintf([files{index} '\nFixed km: (%d) '], index)
    load(['estimates_fixed2/' files{index}])
    loglik_fixed(index, :) = estimator.loglik(nrep);
    fprintf('%.3f\n', estimator.system.fixed.values(end));
    
    fprintf([files{index} '\nVariable km: (%d) '], index)
    load(['estimates_final3/' files{index}])
    loglik_final(index, :) = estimator.loglik(nrep);
    fprintf('%.3f\n\n', estimator.results_GMGTS.b_est(end));
    
    fprintf([files{index} '\nVariable km (GTS): (%d) '], index)
%     load(['estimates_final3/' files{index}])
    loglik_final_GTS(index, :) = estimator.loglik(nrep, "GTS");
    fprintf('%.3f\n\n', estimator.results_GMGTS.b_est(end));
end

save('likelihoods2.mat');

["Marker" "Fixed km" "(SD)" "Flexible km" "(SD)" "Flexible (GTS)" "(SD";
 strrep(strrep(string(files), "EL222_AQTrip-", ""), "EL222_AQTrip_", "") ...
 mean(loglik_fixed')' std(loglik_fixed')' mean(loglik_final')' std(loglik_final')' ...
                                  mean(loglik_final_GTS')' std(loglik_final_GTS')']

%%
plot(estimator, 'States', [1:system.K], 'MaxCells', 100)


% EL222_AQTrip-sfGFP_microfluidic_CFP_1x1s_rep1.mat
% Fixed km: (10) 1 2 3 4 5 6 7 8 9 10 
%  Mean: -177757 SD: 519.99
% 0.099
% EL222_AQTrip-sfGFP_microfluidic_CFP_1x1s_rep1.mat
% Variable km: (10) 1 2 3 4 5 6 7 8 9 10 
%  Mean: -174169 SD: 3969.86
% 0.099
% 
% EL222_AQTrip-sfGFP_microfluidic_CFP_1x1s_rep1.mat
% Variable km (GTS): (10) 1 2 3 4 5 6 7 8 9 10 
%  Mean: 69 SD: 3.17
% 0.099


means = zeros(1, nfiles);
for index = 1:nfiles
    load(['estimates_final3/' files{index}])
    means(index) = estimator.results_GMGTS.b_est(end);
end
[means, order] = sort(means);

names = string(files(order));
for index = 1:nfiles
    names(index) = regexprep(names(index), "EL222_AQTrip[-_]", "");
    names(index) = regexprep(names(index), "_.+", "");
end


ll_fixed = round(mean(loglik_fixed(order, :), 2));
ll_variable = round(mean(loglik_final(order, :), 2));
fprintf('%s\n', strcat(names, repmat(" & ", 12, 1), string(ll_fixed), repmat(" & ", 12, 1), string(ll_variable), repmat(" \\", 12, 1)))





