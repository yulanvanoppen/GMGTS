%% Setup system ------------------------------------------------------------                
clearvars
close all

% system = System('model_repressilator_full.txt', FixedParameters=["DNAT" "kf" "Kd" "m1" "p1"]);
% save('system_repressilator.mat', 'system')
load('system_repressilator.mat')

methods = [];
methods = [methods "GMGTS"];
methods = [methods "GTS"];

seeds = 1;
nseeds = length(seeds);

[ws_GMGTS, ws_GTS, ws_truth, times_GMGTS, times_GTS] = deal(zeros(nseeds, 2));

correlation = eye(6) + diag([.5 0 .5 0 .5], 1) + diag([.5 0 .5 0 .5], -1);
D = .01 * system.k0 .* system.k0' .* correlation;

for first_obs = 1
    generator = Generator(system, N=200, t=0:5:100, error_std=.05, D=D, observed=[1:3 3*first_obs+1:9]);

    for seed = seeds
        rng(seed);
        [data, ground_truth] = generator.generate();
    
        estimator = Estimator(system, data, Methods=methods, Knots=10:10:90);
    
        [GMGTS, GTS] = estimator.estimate();
                                   
        ws_GMGTS(seed, first_obs) = wsdist(GMGTS.b, GMGTS.D, ground_truth.b, ground_truth.D);
        ws_GTS(seed, first_obs) = wsdist(GTS.b, GTS.D, ground_truth.b, ground_truth.D);
        ws_truth(seed, first_obs) = wsdist(mean(ground_truth.beta),cov(ground_truth.beta), ...
                                           ground_truth.b, ground_truth.D);
                                   
        times_GMGTS(seed, first_obs) = sum(estimator.results_GMGTS.time);
        times_GTS(seed, first_obs) = sum(estimator.results_GTS.time);
    end
end




%% Plot --------------------------------------------------------------------
close all


plot(estimator, True=ground_truth, States=1:9, MaxCells=10)

% mkdir('simulation')
% save('simulation/repressilator.mat', 'ws_GMGTS', 'ws_GTS', 'ws_truth', 'times_GMGTS', 'times_GTS');


