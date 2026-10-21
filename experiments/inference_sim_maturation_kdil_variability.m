%% Simulation study: effect of unmodeled dilution rate variability
%
% Generates data from the two-step FP maturation model
% (models/model_maturation_twostep.txt) with a cell-to-cell variable
% dilution rate kdil, and infers the maturation rate distribution with
% GMGTS while treating kdil as fixed (only the final fluorescent state is
% observed). A control experiment uses data generated with fixed kdil.
% The script prints the median and interquartile range of the estimated
% coefficient of variation of km (control first) as a LaTeX table row.
%
% Output: results/simulation/variability_kdil.mat
% Runtime: moderate (20 estimations with N = 500 and N = 100 cells).
%
% Run setup_GMGTS once per MATLAB session before this script.


%% Setup system ------------------------------------------------------------
clearvars
close all

results_dir = fullfile(gmgts_root, 'results', 'simulation');
if ~isfolder(results_dir), mkdir(results_dir); end
wsfile = fullfile(results_dir, 'variability_kdil.mat');

% system = System('model_maturation_onestep.txt', FixedParameters=["kr" "kdr" "kdil" "d"]);
% system_variable = System('model_maturation_onestep.txt', FixedParameters=["kr" "kdr" "d"]);

system = System('model_maturation_twostep.txt', FixedParameters=["kr" "kdr" "kdil" "d"]);
system_variable = System('model_maturation_twostep.txt', FixedParameters=["kr" "kdr" "d"]);

% km_nominal = .025;
km_nominal = .05;

[system.k0(end), system_variable.k0(end)] = deal(km_nominal);


%% Generate data with variable kdil and estimate with fixed kdil -----------
generator = Generator(system_variable, N=500, t=(0:5:100)*(system.K-1), error_std=.01, ...
                      D=diag(([.05 .25 .25] .* system_variable.k0').^2), observed=system.K);            
seeds = 1:10;
nseeds = length(seeds);

km_m = zeros(1, nseeds);                                                    % estimated mean of km
km_lsd = zeros(1, nseeds);                                                  % estimated standard deviation of km

for seed = seeds
    rng(seed);
    [data, ground_truth] = generator.generate();
    
    selection = [2 3];                                                      % ground truth without kdil (for plotting)
    ground_truth.b = ground_truth.b(selection);
    ground_truth.D = ground_truth.D(selection, selection);
    ground_truth.beta = ground_truth.beta(:, selection);

    estimator = Estimator(system, data, Knots=[10 20 60 120], LB=[.001 .001], UB=[10 1]);
    estimator.estimate();

    km_m(seed) = estimator.results_GMGTS.b_est(end);
    km_lsd(seed) = sqrt(estimator.results_GMGTS.D_est(end));
end

plot(estimator, True=ground_truth, MaxCells=10)
 


%% Control: data generated with fixed kdil ---------------------------------
generator = Generator(system, N=100, t=(0:5:100)*(system.K-1), error_std=.01, ...
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


%% LaTeX table row: median [IQR] of the estimated CV of km -----------------
% e.g. 0.203 & \textit{0.169} & \textit{0.209} && 0.202 & \textit{0.189} & \textit{0.221}
quants = quantile([km_sd_fixed' ./ km_m_fixed', km_lsd' ./ km_m'], [.5 .25 .75], 1);
strings = string(num2str(quants(:), '%.3f'));

strcat(strings(1), " & \textit{", strings(2), "} & \textit{", strings(3), "} && ", ... 
       strings(4), " & \textit{", strings(5), "} & \textit{", strings(6), "}")

save(wsfile);
