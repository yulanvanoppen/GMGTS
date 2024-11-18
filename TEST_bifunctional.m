%% Setup system ------------------------------------------------------------
                          
clearvars
close all

model =  'model_bifunctional.txt';
system = ODEIQM(model , 'FixedParameters', ["k1" "k2" "k4" "k6"]);

save('system_bifunctional_measurable.mat', 'system')

load('system_bifunctional_measurable.mat')

first_obs = 5;
dt = 10;
noise_level = .05;
seed = 1;

generator = Generator(system ...                                            % generator setup
                      , 'N', 1000 ...                                        % number of cells
                      , 't', unique([0 2.5 5 dt:dt:100]) ...                % time grid
                      , 'error_std', noise_level ...                        % std of lognormal multiplicative errors
                      , 'D_mult', .25 ...                                   % variance scale
                      , 'observed', first_obs:system.K ...                                 % observed states labels (or indices)
                      );

rng(seed);
[data, ground_truth] = generator.generate();
% plot(generator)


%% Estimate ----------------------------------------------------------------

methods = [];
methods = [methods "GMGTS"];
% methods = [methods "GTS"];

estimator = Estimator(system, data ...                                      % estimator setup
                      , 'Stages', 2 ...                                     % 0: smoothing only, 1: first stage only
                      , 'Methods', methods ...                              % GMGT, GTS, or both
                      , 'MaxIterationsFS', 5 ...
                      , 'ConvergenceTolFS', 1e-12 ...
                      , 'Knots', [0 2.5 40] ...
                  );

estimator.estimate();

% close all
plot(estimator, 'True', ground_truth ...
     , 'States', 1:6 ...
     , 'MaxCells', 7)

