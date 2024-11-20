%% Setup system ------------------------------------------------------------
                          
clearvars
close all

% model =  'model_repressilator_full.txt';
% system = ODEIQM(model, 'FixedParameters', ["DNAT" "kf" "Kd" "m1" "p1"]);
% 
% save('system_repressilator.mat', 'system')

load('system_repressilator.mat')

observed = [1:3 7:9];
% observed = 1:9;
seed = 1;

correlation = eye(6) + diag([.5 0 .5 0 .5], 1) + diag([.5 0 .5 0 .5], -1);
D = .01 * system.k0 .* system.k0' .* correlation;

generator = Generator(system, N=100, t=0:5:100, error_std=.05, D=D, observed=observed);

rng(seed);
[data, ground_truth] = generator.generate();
% plot(generator)


%% Estimate ----------------------------------------------------------------

methods = [];
methods = [methods "GMGTS"];
% methods = [methods "GTS"];

estimator = Estimator(system, data, Stages=2, Methods=methods,...
                      Knots=linspace(0, 100, 14));
rng(seed);
estimator.estimate();

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
% plot(estimator, True=ground_truth, States=1:9, MaxCells=7)



