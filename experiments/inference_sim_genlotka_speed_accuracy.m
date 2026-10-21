%% Simulation study: scalability for generalized Lotka-Volterra models
%
% Compares the accuracy (normalized Wasserstein distance) and computation
% times of GMGTS and GTS on simulated data from 16-state generalized
% Lotka-Volterra models with P = 1, 2, 4, 8 random effects
% (models/model_generalizedLV<P>.txt), for full observation and with the
% first four states hidden, using 10 random seeds per setting.
%
% Output: results/simulation/genLV_speed.mat
% Runtime: long (4 x 10 x 2 = 80 estimations with N = 100 cells).
%
% Run setup_GMGTS once per MATLAB session before this script.


%% Setup system ------------------------------------------------------------
                          
clearvars
close all

results_dir = fullfile(gmgts_root, 'results', 'simulation');
if ~isfolder(results_dir), mkdir(results_dir); end
wsfile = fullfile(results_dir, 'genLV_speed.mat');

P_values = [1 2 4 8];
seeds = 1:10;

[ws_GMGTS, ws_GTS, ws_truth, times_GMGTS, times_GTS] = deal(zeros(length(seeds), length(P_values), 2));

for first_obs = [1 5]                                                       % 1: all states observed, 5: first four states hidden
    for P_idx = 1:length(P_values)
        for seed = seeds
            P = P_values(P_idx);
            P_seed = [P seed]
            
            system = System(sprintf('model_generalizedLV%d.txt', P), FixedParameters=strcat('r', string(1:16)));

            generator = Generator(system, N=100, t=0:20, error_std=.05, D_mult=.25, observed=first_obs:system.K);

            first_idx = min(first_obs, 2);
            rng(100*first_idx + 10*P_idx + seed - 110 - (min(seeds) > 0));
            [data, ground_truth] = generator.generate();


            %% Estimate ----------------------------------------------------

            methods = [];
            methods = [methods "GMGTS"];
            methods = [methods "GTS"];

            estimator = Estimator(system, data, Methods=methods, Knots=[5 10 15]);
            [GMGTS, GTS] = estimator.estimate(true);

            GMGTS_wasserstein = wsdist(GMGTS.b, GMGTS.D, ground_truth.b, ground_truth.D);
            GTS_wasserstein = wsdist(GTS.b, GTS.D, ground_truth.b, ground_truth.D);
            truth_wasserstein = wsdist(mean(ground_truth.beta), cov(ground_truth.beta), ground_truth.b, ground_truth.D);
            
            ws_GMGTS(max(1, seed), P_idx, first_idx) = GMGTS_wasserstein;
            ws_GTS(max(1, seed), P_idx, first_idx) = GTS_wasserstein;
            ws_truth(max(1, seed), P_idx, first_idx) = truth_wasserstein;

            times_GMGTS(max(1, seed), P_idx, first_idx) = sum(estimator.results_GMGTS.time);
            times_GTS(max(1, seed), P_idx, first_idx) = sum(estimator.results_GTS.time);
            
            ws_GMGTS_GTS = [GMGTS_wasserstein GTS_wasserstein]
            times_GMGTS_GTS = [sum(estimator.results_GMGTS.time) sum(estimator.results_GTS.time)]
        end
    end
end

save(wsfile);
