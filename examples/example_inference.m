%% GMGTS example: simulate data from a mixed-effects model and infer it
%
% Generates single-cell time series from the generalized Lotka-Volterra
% model with P = 4 random effects (model_generalizedLV4.txt), of which the
% first four states are hidden, and recovers the random effects
% distribution with GMGTS (optionally also with GTS for comparison).
%
% Requirements: MATLAB R2021a or later, IQM Tools, and the toolboxes listed
% in README.md. Run setup_GMGTS once per MATLAB session before this script.
% Runtime: about a minute (excluding time spent in the interactive app).


%% Setup system ------------------------------------------------------------

P = 4;                                                                      % number of random effects (1, 2, 4, or 8)
model_file = sprintf('model_generalizedLV%d.txt', P);                      % (located in models/)
fixed_param = strcat('r', string(1:16));                                    % growth rates r1,...,r16 are fixed


%% Generate data using the Generator class

system = System(model_file, FixedParameters=fixed_param);                   % create ODE System object

                                                                            % setup Generator for 25 cells at time points 0,1,...,20
generator = Generator(system, N=25, t=0:20, error_std=.1, ...               % with multiplicative noise std 0.1, random effect coefficient 
                      D_mult=.25, observed=5:system.K);                     % of variation 0.25, and the first four states hidden
                      

rng(0);                                                                     % fix random seed for reproducibility
[data, ground_truth] = generator.generate()                                 % generate data and save ground truth
plot(generator)


%% Estimate using the Estimator class

methods = "GMGTS";
% methods = [methods "GTS"];                                                % uncomment to add GTS inference

                                                                            % setup Estimator for the chosen methods
                                                                            % with smoothing using the interactive app
estimator = Estimator(system, data, Methods=methods, InteractiveSmoothing=true);

                                                                            % in case GTS is also used, optionally replace with
estimates = estimator.estimate()                                            % [GMGTSest, GTSest] = estimator.estimate()

                                                                            % plot estimates/predictions versus ground truth
plot(estimator, True=ground_truth, States=9:system.K, MaxCells=7)           % for the last 8 states and the first 7 cells


%% Equivalent: estimate using the GMGTS() function
                                                                            % utility function that takes care of System/Estimator instantiation
estimates = GMGTS(model_file, data, FixedParameters=fixed_param, Methods=methods, InteractiveSmoothing=true)


%% Equivalent: estimate using the GMGTS() function without a data struct
y = data.traces;
t = data.t;
observed = data.observed;                                                   % equivalent with separate measurement arguments
init = data.init;                                                           % omitting optional positional arguments t, observed, init
                                                                            % leads to defaults being used instead
estimates = GMGTS(model_file, y, t, observed, init, FixedParameters=fixed_param, ...
                  Methods=methods, InteractiveSmoothing=true, Plot=false)   % disable plots
