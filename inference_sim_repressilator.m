%% Setup system ------------------------------------------------------------                
clearvars
close all

name = 'repressilator';
model =  'model_repressilator_full.txt';
system = System(name, model, 'FixedParameters', ["DNAT" "kf" "Kd" "m1" "p1"]);

save('system_repressilator.mat', 'system')

load('system_repressilator.mat')


methods = [];
methods = [methods "GMGTS"];
methods = [methods "GTS"];
initial = system.k0';

                  
seeds = 1:10;
nseeds = length(seeds);

GMGTS_distances_ws = zeros(nseeds, 2);
GTS_distances_ws = zeros(nseeds, 2);
truth_distances_ws = zeros(nseeds, 2);

times_GMGTS = zeros(nseeds, 2);
times_GTS = zeros(nseeds, 2);

datas = cell(nseeds, 2);
truths = cell(nseeds, 2);
estimators = cell(nseeds, 2);

correlation = eye(6) + diag([.5 0 .5 0 .5], 1) + diag([.5 0 .5 0 .5], -1);
D = .01 * system.k0 .* system.k0' .* correlation;

for first_obs = [1 2]
    generator = Generator(system ...                                            % generator setup
                      , 'N', 100 ...                                        % number of cells
                      , 't', [0:5:100] ...                                  % time grid
                      , 'error_std', .05 ...                               % std of lognormal multiplicative errors
                      , 'D', D ...                                   % variance scale
                      , 'observed', [1:3 3*first_obs+1:9] ...             % observed states labels (or indices)
                      , 'varying', 1:system.P ...                           % variable parameter labels (or indices)
                      );

    for seed = seeds
        rng(seed);
        seed = max(1, seed);
        generator.generate();
    
        generated = generator.data;
        [datas{seed, first_obs}, truths{seed, first_obs}] = obfuscate(generated);
    
        datas{seed, first_obs}.beta = truths{seed, first_obs}.beta;
        datas{seed, first_obs}.original = truths{seed, first_obs}.original;
        datas{seed, first_obs}.doriginal = truths{seed, first_obs}.doriginal;
    
    
        estimators{seed, first_obs} = Estimator(datas{seed, first_obs}, system ...
                                     , 'Stages', 2 ...
                                     , 'Methods', methods ...
                                     );
    
        estimators{seed, first_obs}.estimate();
                                   
        GMGTS_distances_ws(seed, first_obs) = wsdist(estimators{seed, first_obs}.results_GMGTS.b_est, estimators{seed, first_obs}.results_GMGTS.D_GTS, ...
                                       truths{seed, first_obs}.b, truths{seed, first_obs}.D);
        GTS_distances_ws(seed, first_obs) = wsdist(estimators{seed, first_obs}.results_GTS.b_est, estimators{seed, first_obs}.results_GTS.D_GTS, ...
                                     truths{seed, first_obs}.b, truths{seed, first_obs}.D);
        truth_distances_ws(seed, first_obs) = wsdist(mean(truths{seed, first_obs}.beta),cov(truths{seed, first_obs}.beta), ...
                                       truths{seed, first_obs}.b, truths{seed, first_obs}.D);
                                   
        times_GMGTS(seed, first_obs) = sum(estimators{seed, first_obs}.results_GMGTS.time);
        times_GTS(seed, first_obs) = sum(estimators{seed, first_obs}.results_GTS.time);
                                   
        fprintf('GMGTS Wasserstein distance: %.4f\n', GMGTS_distances_ws(seed, first_obs));
        
    end
end




%% Plot --------------------------------------------------------------------
close all

obs_idx = 2;

[~, median_idx_partial_GTS] = min(abs(GTS_distances_ws(:, obs_idx) - median(GTS_distances_ws(:, obs_idx))))
[~, median_idx_partial] = min(abs(GMGTS_distances_ws(:, obs_idx) - median(GMGTS_distances_ws(:, obs_idx))))
median_dist_partial = GMGTS_distances_ws(median_idx_partial, obs_idx)

plot(estimators{median_idx_partial, obs_idx}, 'True', truths{median_idx_partial, obs_idx} ...
     , 'States', [1:4] ...
     , 'MaxCells', 10)
 

mkdir('simulation')
wsfile = 'simulation/repressilator.mat';
save(wsfile);


