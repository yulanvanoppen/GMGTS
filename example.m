%% Setup system ------------------------------------------------------------
                          
clearvars
close all

P = 8;
model_file = sprintf('model_generalizedLV%d.txt', P);



%% Generate data using the Generator class

system = System(model_file, FixedParameters=strcat('r', string(1:16)));     % create ODE System object

                                                                            % setup Generator for 100 cells at time points 0,1,...,20
generator = Generator(system, N=100, t=0:20, error_std=.05, ...             % with multiplicative noise std 0.05, random effect coefficient 
                      D_mult=.25, observed = 5:system.K);                   % of variation 0.25, and the first four states hidden
                      

rng(0);
[data, ground_truth] = generator.generate();                                % generate data and save ground truth
plot(generator)


%% Estimate using the Estimator class

methods = "GMGTS";
% methods = [methods "GTS"];                                                % uncomment to add GTS inference

                                                                            % setup Estimator for the chosen methods
                                                                            % with smoothing using the interactive app
estimator = Estimator(model_file, data.traces, data.t, Methods=methods, InteractiveSmoothing=true, FixedParameters=strcat('r', string(1:16)));
estimator.estimate()                                                     
                                                                            % plot estimates/predictions versus ground truth
plot(estimator, True=ground_truth, States=9:system.K, MaxCells=7)           % for the last 8 states and the first 7 cells

