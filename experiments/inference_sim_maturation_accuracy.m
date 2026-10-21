%% Simulation study: accuracy and speed for the FP maturation model
%
% Compares GMGTS and GTS on simulated data from the one-step fluorescent
% protein maturation model (models/model_maturation_onestep.txt) for
% varying sampling intervals (dt), multiplicative noise levels, and
% observation scenarios (both states observed vs. only the fluorescent
% state observed), using 10 random seeds per setting. Accuracy is measured
% by the (normalized) Wasserstein and Hellinger distances between the
% estimated and the data-generating random effects distributions.
%
% Output: results/simulation/maturation_accuracy.mat (saved after every
% completed dt value, so partial results are kept).
% Runtime: long (4 x 6 x 10 x 2 = 480 estimations with N = 500 cells).
%
% Run setup_GMGTS once per MATLAB session before this script.


%% Setup system ------------------------------------------------------------

clearvars
close all

results_dir = fullfile(gmgts_root, 'results', 'simulation');
if ~isfolder(results_dir), mkdir(results_dir); end
wsfile = fullfile(results_dir, 'maturation_accuracy.mat');

system = System('model_maturation_onestep.txt', FixedParameters=["kr" "kdr" "kdil" "d"]);

dt_values = [5 10 20 40];
noise_levels = [.001 .005 .01 .02 .05 .1];
seeds = 1:10;

[ws_GMGTS, ws_GTS, ws_truth, wsu_GMGTS, wsu_GTS, wsu_truth, ...
 he_GMGTS, he_GTS, he_truth, times_GMGTS, times_GTS] = deal(zeros(length(seeds), length(noise_levels), length(dt_values), 2));

for first_obs = [1 2]                                                       % 1: both states observed, 2: only F observed
    for dt_idx = 1:length(dt_values) 
        for noise_idx = 1:length(noise_levels)
            for seed = seeds
                dt = dt_values(dt_idx);
                noise_level = noise_levels(noise_idx);
                dt_noise_seed = [dt noise_level seed]
                
                generator = Generator(system, N=500, t=unique([0 5 10 20 dt:dt:200]), error_std=noise_level, ...
                                      D_mult=.25, observed=first_obs:system.K);

                rng(1000*first_obs + 100*dt_idx + 10*noise_idx + seed - 1110 - (min(seeds) > 0));
                
                [data, ground_truth] = generator.generate();
                

                %% Estimate ------------------------------------------------

                methods = [];
                methods = [methods "GMGTS"];
                methods = [methods "GTS"]; 

                estimator = Estimator(system, data, Methods=methods, Knots=[10 20 60 120], LB=[.001 .001], UB=[20 1], MaxIterationsFS=10);
                [GMGTS, GTS] = estimator.estimate(true);

                GMGTS_wasserstein = wsdist(GMGTS.b, GMGTS.D, ground_truth.b, ground_truth.D);
                GMGTS_wasserstein_unscaled = wsdist(GMGTS.b, GMGTS.D, ground_truth.b, ground_truth.D, false);
                GMGTS_hellinger = hedist(GMGTS.b, GMGTS.D, ground_truth.b, ground_truth.D);
                ws_GMGTS(max(1, seed), noise_idx, dt_idx, first_obs) = GMGTS_wasserstein;
                wsu_GMGTS(max(1, seed), noise_idx, dt_idx, first_obs) = GMGTS_wasserstein_unscaled;
                he_GMGTS(max(1, seed), noise_idx, dt_idx, first_obs) = GMGTS_hellinger;
                
                truth_wasserstein = wsdist(mean(ground_truth.beta),cov(ground_truth.beta), ground_truth.b, ground_truth.D);
                truth_wasserstein_unscaled = wsdist(mean(ground_truth.beta),cov(ground_truth.beta), ground_truth.b, ground_truth.D, false);
                truth_hellinger = hedist(mean(ground_truth.beta),cov(ground_truth.beta), ground_truth.b, ground_truth.D);
                ws_truth(max(1, seed), noise_idx, dt_idx, first_obs) = truth_wasserstein;
                wsu_truth(max(1, seed), noise_idx, dt_idx, first_obs) = truth_wasserstein_unscaled;
                he_truth(max(1, seed), noise_idx, dt_idx, first_obs) = truth_hellinger;
                
                GTS_wasserstein = wsdist(GTS.b, GTS.D, ground_truth.b, ground_truth.D);
                GTS_wasserstein_unscaled = wsdist(GTS.b, GTS.D, ground_truth.b, ground_truth.D, false);
                GTS_hellinger = hedist(GTS.b, GTS.D, ground_truth.b, ground_truth.D);
                ws_GTS(max(1, seed), noise_idx, dt_idx, first_obs) = GTS_wasserstein;
                wsu_GTS(max(1, seed), noise_idx, dt_idx, first_obs) = GTS_wasserstein_unscaled;
                he_GTS(max(1, seed), noise_idx, dt_idx, first_obs) = GTS_hellinger;
                
                ws_GMGTS_GTS = [GMGTS_wasserstein GTS_wasserstein]
                times_GMGTS_GTS = [sum(estimator.results_GMGTS.time) sum(estimator.results_GTS.time)]
                
                times_GMGTS(max(1, seed), noise_idx, dt_idx, first_obs) = times_GMGTS_GTS(1);
                times_GTS(max(1, seed), noise_idx, dt_idx, first_obs) = times_GMGTS_GTS(2);
            end
        end
        save(wsfile);                                                       % checkpoint
    end
end

save(wsfile);
