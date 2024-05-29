%% Setup system ------------------------------------------------------------
clearvars
close all


name = 'maturation_fluorescence'; 
model =  'model_maturation_onestep.txt';
% model =  'model_maturation_twostep.txt';
system = ODEIQM(name, model, 'FixedParameters', ["kr" "kdr" "kdil" "d"]);

save('system_maturation.mat', 'system')


name = 'maturation_fluorescence'; 
model =  'model_maturation_onestep.txt';
% model =  'model_maturation_twostep.txt';
system = ODEIQM(name, model, 'FixedParameters', ["kr" "kdr" "d"]);

save('system_maturation_variable.mat', 'system')




% Generate data -----------------------------------------------------------
load('system_maturation_variable.mat')
generator = Generator(system ...                                            % generator setup
                      , 'N', 100 ...                                        % number of cells
                      , 't', [0:5:100] ...                                  % time grid
                      , 'error_std', .01 ...
                      , 'D', diag(([.05 .25 .25] .* system.k0').^2) ...
                      , 'observed', system.K ...                            % observed states labels (or indices)
                      , 'varying', 1:system.P ...                           % variable parameter labels (or indices)
                      );
load('system_maturation.mat')      

methods = [];
methods = [methods "GMGTS"];
                  
seeds = 1:10;
nseeds = length(seeds);

GMGTS_distances = zeros(1, nseeds);
GTS_distances = zeros(1, nseeds);
truth_distances = zeros(1, nseeds);

km_m = zeros(1, nseeds);
km_sd = zeros(1, nseeds);

datas = cell(1, nseeds);
truths = cell(1, nseeds);
estimators = cell(1, nseeds);


for seed = seeds
    rng(seed);
    seed = max(1, seed);

    generator.generate();

    generated = generator.data;
    [datas{seed}, truths{seed}] = obfuscate(generated);
    
    selection = [2 3];
    datas{seed}.varying = [1 2];
    truths{seed}.b = truths{seed}.b(selection);
    truths{seed}.D = truths{seed}.D(selection, selection);
    truths{seed}.beta = truths{seed}.beta(:, selection);


    estimators{seed} = Estimator(datas{seed}, system ...
                                 , 'Stages', 2 ...
                                 , 'Methods', methods ...
                                 , 'Knots', [12.5 25 50] ...
                                 , 'LB', [.001 .001] ...
                                 , 'UB', [10 1] ...
                                 );

    estimators{seed}.estimate();
    
    GMGTS_distances(seed) = wsdist(estimators{seed}.results_GMGTS.b_est, estimators{seed}.results_GMGTS.D_GTS, ...
                                   truths{seed}.b, truths{seed}.D);
    truth_distances(seed) = wsdist(mean(truths{seed}.beta),cov(truths{seed}.beta), ...
                                   truths{seed}.b, truths{seed}.D);

    km_m(seed) = estimators{seed}.results_GMGTS.b_est(end);
    km_sd(seed) = sqrt(estimators{seed}.results_GMGTS.D_GTS(end, end));
end



[~, median_idx] = min(abs(GMGTS_distances - median(GMGTS_distances)))
median_dist = GMGTS_distances(median_idx)

plot(estimators{median_idx}, 'True', truths{median_idx} ...
     , 'States', [1:system.K] ...
     , 'MaxCells', 20)
 
km_sd' ./ km_m'
 


%% Control
load('system_maturation_delay.mat')
generator = Generator(system ...                                            % generator setup
                      , 'N', 100 ...                                        % number of cells
                      , 't', [0:10:200] ...                                  % time grid
                      , 'error_std', .01 ...
                      , 'D', diag(([.25 .25] .* system.k0').^2) ...
                      , 'observed', system.K ...                            % observed states labels (or indices)
                      , 'varying', 1:system.P ...                           % variable parameter labels (or indices)
                      );

GMGTS_distances_fixed = zeros(1, nseeds);
GTS_distances_fixed = zeros(1, nseeds);
truth_distances_fixed = zeros(1, nseeds);

km_m_fixed = zeros(1, nseeds);
km_sd_fixed = zeros(1, nseeds);

datas_fixed = cell(1, nseeds);
truths_fixed = cell(1, nseeds);
estimators_fixed = cell(1, nseeds);

for seed = seeds
    rng(seed);
    seed = max(1, seed);

    generator.generate();

    generated = generator.data;
    [datas_fixed{seed}, truths_fixed{seed}] = obfuscate(generated);
   

    estimators_fixed{seed} = Estimator(datas_fixed{seed}, system ...
                                 , 'Stages', 2 ...
                                 , 'Methods', methods ...
                                 , 'Knots', [12.5 25 50] ...
                                 , 'LB', [.001 .001] ...
                                 , 'UB', [10 1] ...
                                 );

    estimators_fixed{seed}.estimate();
    
    GMGTS_distances_fixed(seed) = wsdist(estimators_fixed{seed}.results_GMGTS.b_est, estimators_fixed{seed}.results_GMGTS.D_GTS, ...
                                   truths_fixed{seed}.b, truths_fixed{seed}.D);
    truth_distances(seed) = wsdist(mean(truths_fixed{seed}.beta),cov(truths_fixed{seed}.beta), ...
                                   truths_fixed{seed}.b, truths_fixed{seed}.D);

    km_m_fixed(seed) = estimators_fixed{seed}.results_GMGTS.b_est(end);
    km_sd_fixed(seed) = sqrt(estimators_fixed{seed}.results_GMGTS.D_GTS(end, end));
end



[~, median_idx] = min(abs(GMGTS_distances_fixed - median(GMGTS_distances_fixed)))
median_dist = GMGTS_distances_fixed(median_idx)

plot(estimators_fixed{median_idx}, 'True', truths_fixed{median_idx} ...
     , 'States', [1:system.K] ...
     , 'MaxCells', 20)

[km_sd_fixed' ./ km_m_fixed', km_sd' ./ km_m']
quantile([km_sd_fixed' ./ km_m_fixed', km_sd' ./ km_m'], [.25 .5 .75], 1)



%% LaTeX
mkdir('simulation')
wsfile = 'simulation/variability_kdil.mat';
save(wsfile);

