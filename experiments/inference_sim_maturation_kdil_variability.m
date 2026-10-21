%% Setup system ------------------------------------------------------------
clearvars
close all

% system = System('model_maturation_onestep.txt', 'FixedParameters', ["kr" "kdr" "kdil" "d"]);
% system_variable = System('model_maturation_onestep.txt', 'FixedParameters', ["kr" "kdr" "d"]);
% save('system_maturation_variable.mat', 'system', 'system_variable')

system = System('model_maturation_twostep.txt', 'FixedParameters', ["kr" "kdr" "kdil" "d"]);
system_variable = System('model_maturation_twostep.txt', 'FixedParameters', ["kr" "kdr" "d"]);
save('system_maturation_variable.mat', 'system', 'system_variable')

load('system_maturation_variable.mat')

% km_nominal = .025;
km_nominal = .05;

[system.k0(end), system_variable.k0(end)] = deal(km_nominal);


% Generate data -----------------------------------------------------------
generator = Generator(system_variable, N=500, t=(0:5:100)*(system.K-1) ,error_std=.01, ...
                      D=diag(([.05 .25 .25] .* system_variable.k0').^2), observed=system.K);            
seeds = 1:10;
nseeds = length(seeds);

km_m = zeros(1, nseeds);
km_sd = zeros(1, nseeds);

for seed = seeds
    rng(seed);
    [data, ground_truth] = generator.generate();
    
    selection = [2 3];
    ground_truth.b = ground_truth.b(selection);
    ground_truth.D = ground_truth.D(selection, selection);
    ground_truth.beta = ground_truth.beta(:, selection);

    estimator = Estimator(system, data, Knots=[10 20 60 120], LB=[.001 .001], UB=[10 1]);
    estimator.estimate();

    km_m(seed) = estimator.results_GMGTS.b_est(end);
    km_lsd(seed) = sqrt(estimator.results_GMGTS.D_est(end));
end

plot(estimator, True=ground_truth, MaxCells=10)
 


%% Control
generator = Generator(system, N=100, t=(0:5:100)*(system.K-1) ,error_std=.01, ...
                      D=diag(([.25 .25] .* system.k0').^2), observed=system.K);

km_m_fixed = zeros(1, nseeds);
km_sd_fixed = zeros(1, nseeds);

for seed = seeds
    rng(seed);
    [data, ground_truth] = generator.generate();

    estimator = Estimator(system, data, Knots=[10 20 60 120], LB=[.001 .001], UB=[10 1]);
    estimator.estimate();

    km_m_fixed(seed) = estimator.results_GMGTS.b_est(end);
    km_sd_fixed(seed) = sqrt(estimator.results_GMGTS.D_est(end));
end

% & 0.203 & \textit{0.169} & \textit{0.209} && 0.202 & \textit{0.189} & \textit{0.221}
quants = quantile([km_sd_fixed' ./ km_m_fixed', km_lsd' ./ km_m'], [.5 .25 .75], 1);
strings = string(num2str(quants(:), '%.3f'));

strcat(strings(1), " & \textit{", strings(2), "} & \textit{", strings(3), "} && ", ... 
       strings(4), " & \textit{", strings(5), "} & \textit{", strings(6), "}")



%% LaTeX
mkdir('simulation')
wsfile = 'simulation/variability_kdil.mat';
save(wsfile);

