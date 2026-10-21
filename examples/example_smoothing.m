%% GMGTS example: interactive B-spline smoothing only
%
% Loads a precomputed System object, measurements, and ground truth from
% example_data.mat and opens the interactive smoothing app (Stages=0, so
% only the smoothing step of GMGTS is carried out). IQM Tools is not
% required for this example.
%
% Run setup_GMGTS once per MATLAB session before this script.


%% Load example system and data

load(fullfile(gmgts_root, 'examples', 'example_data.mat'))                 % loads system, data, ground_truth


%% Smooth using the Estimator class                                         % (avoids System construction using IQM Tools)

methods = "GMGTS";
                                                                            % setup Estimator for smoothing only (Stages=0)
                                                                            % using the interactive app
estimator = Estimator(system, data, Methods=methods, Stages=0, InteractiveSmoothing=true);

estimates = estimator.estimate()                                            % returns smoothed states and gradients

                                                                            % plot smoothed measurements versus ground truth
plot(estimator, True=ground_truth, States=1:system.K, MaxCells=7)           % for all states and the first 7 cells
