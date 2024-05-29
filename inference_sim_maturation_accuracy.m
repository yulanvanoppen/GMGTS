%% Setup system ------------------------------------------------------------

clearvars
close all

name = 'maturation_fluorescence'; 
model =  'model_maturation_onestep.txt';
system = ODEIQM(name, model, 'FixedParameters', ["kr" "kdr" "kdil" "d"]);

% save('system_maturation_delay.mat', 'system')

% load('system_maturation_delay.mat')


dt_values = [5 10 20 40];
noise_levels = [.001 .005 .01 .02 .05 .1];
seeds = 1:10;

wasserstein_GMGTS = zeros(length(seeds), length(noise_levels), length(dt_values), 2);
wasserstein_GTS = zeros(length(seeds), length(noise_levels), length(dt_values), 2);
wasserstein_truth = zeros(length(seeds), length(noise_levels), length(dt_values), 2);

times_GMGTS = zeros(length(seeds), length(noise_levels), length(dt_values), 2);
times_GTS= zeros(length(seeds), length(noise_levels), length(dt_values), 2);

for first_obs = 2
    for dt_idx = 1:length(dt_values) 
        for noise_idx = 1:length(noise_levels)
            for seed = seeds
                dt = dt_values(dt_idx)
                noise_level = noise_levels(noise_idx)
                seed
                
                generator = Generator(system ...                                            % generator setup
                                      , 'N', 100 ...                                        % number of cells
                                      , 't', unique([0 5 10 20 dt:dt:200]) ...                                  % time grid
                                      , 'error_std', noise_level ...                               % std of lognormal multiplicative errors
                                      , 'D_mult', .25 ...
                                      , 'observed', first_obs:system.K ...                            % observed states labels (or indices)
                                      , 'varying', 1:system.P ...                           % variable parameter labels (or indices)
                                      );

                rng(1000*first_obs + 100*dt_idx + 10*noise_idx + seed - 1110 - (min(seeds) > 0));                                 % fix random number generator
                generator.generate();                                                       % generate data and package

                generated = generator.data;
                [data, ground_truth] = obfuscate(generated);
                
                
                %% Estimate ----------------------------------------------------------------

                methods = [];
                methods = [methods "GMGTS"];
                methods = [methods "GTS"];

                estimator = Estimator(data, system ...                                      % estimator setup
                                      , 'Stages', 2 ...                                     % 0: smoothing only, 1: first stage only
                                      , 'Methods', methods ...                              % GMGT, GTS, or both
                                      , 'Knots', [10 50] ...
                                      , 'LB', [.001 .01] ...
                                      , 'UB', [1 1] ...
                                      );

                estimator.estimate('silent');                                                       % estimate parameters
                
                close all
                plot(estimator, 'True', ground_truth ...
                     , 'States', 1:2 ...
                     , 'MaxCells', 20)
                [sqrt(diag(cov(estimator.results_GMGTS.beta_fs-ground_truth.beta))) sqrt(diag(mean(estimator.results_GMGTS.variances_beta_fs, 3)))]

                GMGTS_wasserstein = wsdist(estimator.results_GMGTS.b_est, estimator.results_GMGTS.D_GTS, ground_truth.b, ground_truth.D)
                wasserstein_GMGTS(max(1, seed), noise_idx, dt_idx, first_obs) = GMGTS_wasserstein;
                
                truth_wasserstein = wsdist(mean(ground_truth.beta),cov(ground_truth.beta), ground_truth.b, ground_truth.D);
                wasserstein_truth(max(1, seed), noise_idx, dt_idx, first_obs) = truth_wasserstein;
                
                GTS_wasserstein = wsdist(estimator.results_GTS.b_est, estimator.results_GTS.D_GTS, ground_truth.b, ground_truth.D)
                wasserstein_GTS(max(1, seed), noise_idx, dt_idx, first_obs) = GTS_wasserstein;
                
                times_GMGTS_GTS = [sum(estimator.results_GMGTS.time) sum(estimator.results_GTS.time)]
                
                times_GMGTS(max(1, seed), noise_idx, dt_idx, first_obs) = sum(estimator.results_GMGTS.time);
                times_GTS(max(1, seed), noise_idx, dt_idx, first_obs) = sum(estimator.results_GTS.time);
            end
        end
    end
end


mkdir('simulation')
wsfile = 'simulation/maturation_accuracy.mat';
save(wsfile);


