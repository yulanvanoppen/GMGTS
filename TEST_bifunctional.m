%% Setup system ------------------------------------------------------------
                          
clearvars
close all

% system = System('model_bifunctional2.txt', FixedParameters=["k1" "k2" "k4" "k6"]);
% save('system_bifunctional_measurable.mat', 'system')
load('system_bifunctional_measurable.mat')

first_obs = 1;
dt = 2.5;
noise_level = .01;
seed = 2;

generator = Generator(system, N=10, t=unique([0 2.5 5 dt:dt:100]), error_std=noise_level, ...
                      D_mult=.25, observed=first_obs:system.K);
rng(seed);
[data, ground_truth] = generator.generate();
% plot(generator)


%% Estimate ----------------------------------------------------------------

methods = [];
methods = [methods "GMGTS"];
% methods = [methods "GTS"];

estimator = Estimator(system, data, Stages=2, Methods=methods, ...
                      TestConvergence=true);
rng(seed);
out = estimator.estimate();

[GMGTS_hellinger, GMGTS_wasserstein, GTS_hellinger, GTS_wasserstein] = deal(0); %#ok<ASGLU>

GMGTS_est = estimator.results_GMGTS;
GMGTS_hellinger = hedist(GMGTS_est.b_est, GMGTS_est.D_est, ground_truth.b, ground_truth.D);
GMGTS_wasserstein = wsdist(GMGTS_est.b_est, GMGTS_est.D_est, ground_truth.b, ground_truth.D);

% GTS_est = estimator.results_GTS;
% GTS_hellinger = hedist(GTS_est.b_est, GTS_est.D_est, ground_truth.b, ground_truth.D);
% GTS_wasserstein = wsdist(GTS_est.b_est, GTS_est.D_est, ground_truth.b, ground_truth.D);

hellinger_GMGTS_GTS = [GMGTS_hellinger GTS_hellinger]
wasserstein_GMGTS_GTS = [GMGTS_wasserstein GTS_wasserstein]

close all
plot(estimator, True=ground_truth, States=1:6, MaxCells=7)


