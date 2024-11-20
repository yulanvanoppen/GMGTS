%% Setup system ------------------------------------------------------------
                          
clearvars

% name = 'generalizedLV';
% model =  'model_generalizedLV1.txt';
% system = System(name, model, 'FixedParameters', strcat('r', string(1:16)));
% save('system_generalizedLV1.mat', 'system')
% 
% name = 'generalizedLV';
% model =  'model_generalizedLV2.txt';
% system = System(name, model, 'FixedParameters', strcat('r', string(1:16)));
% save('system_generalizedLV2.mat', 'system')
% 
% name = 'generalizedLV';
% model =  'model_generalizedLV4.txt';
% system = System(name, model, 'FixedParameters', strcat('r', string(1:16)));
% save('system_generalizedLV4.mat', 'system')
% 
% name = 'generalizedLV';
% model =  'model_generalizedLV8.txt';
% system = System(name, model, 'FixedParameters', strcat('r', string(1:16)));
% save('system_generalizedLV8.mat', 'system')

% load('system_generalizedLV1.mat')
% load('system_generalizedLV2.mat')
% load('system_generalizedLV4.mat')
% load('system_generalizedLV8.mat')

P_values = [1 2 4 8];
seeds = 1:10;

accuracies_GMGTS = zeros(length(seeds), length(P_values), 2);
accuracies_GTS = zeros(length(seeds), length(P_values), 2);
accuracies_truth = zeros(length(seeds), length(P_values), 2);

times_GMGTS = zeros(length(seeds), length(P_values), 2);
times_GTS= zeros(length(seeds), length(P_values), 2);

for first_obs = 2
    for P_idx = 1:length(P_values)
        for seed = seeds
            P = P_values(P_idx)
            seed
            
            load(sprintf('system_generalizedLV%d.mat', P));

            if first_obs == 1, observed = 1:system.K; else, observed = 5:system.K; end

            generator = Generator(system ...                                            % generator setup
                                  , 'N', 100 ...                                        % number of cells
                                  , 't', [0:1:20] ...                                   % time grid
                                  , 'error_std', .05 ...                                % std of lognormal multiplicative errors
                                  , 'D_mult', .25 ...
                                  , 'observed', observed ...                            % observed states labels (or indices)
                                  , 'varying', 1:system.P ...                           % variable parameter labels (or indices)
                                  );

            rng(100*first_obs + 10*P_idx + seed - 110 - (min(seeds) > 0));              % fix random number generator
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
                                  );

            estimator.estimate('silent');                                               % estimate parameters

            GMGTS_distance = wsdist(estimator.results_GMGTS.b_est, estimator.results_GMGTS.D_GTS, mean(ground_truth.beta), cov(ground_truth.beta))
            truth_distance = wsdist(mean(ground_truth.beta), cov(ground_truth.beta), ground_truth.b, ground_truth.D)

            accuracies_GMGTS(max(1, seed), P_idx, first_obs) = GMGTS_distance;
            accuracies_truth(max(1, seed), P_idx, first_obs) = truth_distance;

            times_GMGTS(max(1, seed), P_idx, first_obs) = sum(estimator.results_GMGTS.time);

            GTS_distance = wsdist(estimator.results_GTS.b_est, estimator.results_GTS.D_GTS, mean(ground_truth.beta), cov(ground_truth.beta))
            accuracies_GTS(max(1, seed), P_idx, first_obs) = GTS_distance;
            times_GTS(max(1, seed), P_idx, first_obs) = sum(estimator.results_GTS.time);
        end
    end
end


mkdir('simulation')
wsfile = 'simulation/genLV_speed.mat';
save(wsfile);


