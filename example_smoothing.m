%% Load example system and data

load('example_data.mat')


%% Estimate using the Estimator class                                       % (avoid System construction using IQM Tools)

methods = "GMGTS";
% methods = [methods "GTS"];                                                % uncomment to add GTS inference

                                                                            % setup Estimator for the chosen methods
                                                                            % with smoothing using the interactive app
estimator = Estimator(system, data, Methods=methods, Stages=0, InteractiveSmoothing=true);

                                                                            % in case GTS is also used, optionally replace with
estimates = estimator.estimate()                                            % [GMGTSest, GTSest] = estimator.estimate()

                                                                            % plot estimates/predictions versus ground truth
plot(estimator, True=ground_truth, States=1:system.K, MaxCells=7)           % for the last 8 states and the first 7 cells