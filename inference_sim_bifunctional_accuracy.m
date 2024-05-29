%% Setup system ------------------------------------------------------------
                          
clearvars
close all

name = 'bifunctional_TCS';
model =  'model_bifunctional.txt';
system = ODEIQM(name, model , 'FixedParameters', ["k1" "k2" "k4" "k6"]);

% save('system_bifunctional_measurable.mat', 'system')

% load('system_bifunctional_measurable.mat')


dt_values = [2.5 5 10 20];
noise_levels = [.001 .005 .01 .02 .05 .1];
seeds = 1:10;

wasserstein_GMGTS = zeros(length(seeds), length(noise_levels), length(dt_values), 2);
wasserstein_GTS = zeros(length(seeds), length(noise_levels), length(dt_values), 2);
wasserstein_truth = zeros(length(seeds), length(noise_levels), length(dt_values), 2);

times_GMGTS = zeros(length(seeds), length(noise_levels), length(dt_values), 2);
times_GTS= zeros(length(seeds), length(noise_levels), length(dt_values), 2);


for first_obs = [1 5]
    for dt_idx = 1:length(dt_values)
        for noise_idx = 1:length(noise_levels)
            for seed = seeds
                dt = dt_values(dt_idx)
                noise_level = noise_levels(noise_idx)
                seed

                generator = Generator(system ...                                            % generator setup
                                      , 'N', 100 ...                                        % number of cells
                                      , 't', unique([0 2.5 5 dt:dt:100]) ...                % time grid
                                      , 'error_std', noise_level ...                        % std of lognormal multiplicative errors
                                      , 'D_mult', .25 ...                                   % variance scale
                                      , 'observed', first_obs:system.K ...                                 % observed states labels (or indices)
                                      , 'varying', 1:system.P ...                           % variable parameter labels (or indices)
                                      );

                first_idx = min(first_obs, 2);
                
                rng(1000*first_idx + 100*dt_idx + 10*noise_idx + seed - 1110 - (min(seeds) > 0));                                 % fix random number generator
                generator.generate();                                                       % generate data and package

                generated = generator.data;
                [data, ground_truth] = obfuscate(generated);
                data.beta = ground_truth.beta;

                %% Estimate ----------------------------------------------------------------

                methods = [];
                methods = [methods "GMGTS"];
                methods = [methods "GTS"];

                estimator = Estimator(data, system ...                                      % estimator setup
                                      , 'Stages', 2 ...                                     % 0: smoothing only, 1: first stage only
                                      , 'Methods', methods ...                              % GMGT, GTS, or both
                                      , 'Knots', [0 2.5 40] ...
                                      , 'Linear', [5 100] ...
                                      );
                
                estimator.estimate('silent');
                
                close all
                plot(estimator, 'True', ground_truth ...
                     , 'States', 1:6 ...
                     , 'MaxCells', 20)
                 [sqrt(diag(cov(estimator.results_GMGTS.beta_fs-ground_truth.beta))) sqrt(diag(mean(estimator.results_GMGTS.variances_beta_fs, 3)))]

                GMGTS_wasserstein = wsdist(estimator.results_GMGTS.b_est, estimator.results_GMGTS.D_GTS, ground_truth.b, ground_truth.D)
                wasserstein_GMGTS(max(1, seed), noise_idx, dt_idx, first_idx) = GMGTS_wasserstein;
                
                truth_wasserstein = wsdist(mean(ground_truth.beta),cov(ground_truth.beta), ground_truth.b, ground_truth.D);
                wasserstein_truth(max(1, seed), noise_idx, dt_idx, first_idx) = truth_wasserstein;
                
                GTS_wasserstein = wsdist(estimator.results_GTS.b_est, estimator.results_GTS.D_GTS, ground_truth.b, ground_truth.D)
                wasserstein_GTS(max(1, seed), noise_idx, dt_idx, first_idx) = GTS_wasserstein;
                
                times_GMGTS_GTS = [sum(estimator.results_GMGTS.time) sum(estimator.results_GTS.time)]
                
                times_GMGTS(max(1, seed), noise_idx, dt_idx, first_idx) = sum(estimator.results_GMGTS.time);
                times_GTS(max(1, seed), noise_idx, dt_idx, first_idx) = sum(estimator.results_GTS.time);
            end
        end
    end
end


mkdir('simulation')
wsfile = 'simulation/bifunctional_accuracy.mat';
save(wsfile);


