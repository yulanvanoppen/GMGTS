%% Setup system ------------------------------------------------------------
                          
clearvars
close all

model =  'model_maturation_onestep.txt';
system = ODEIQM(model, 'FixedParameters', ["kr" "kdr" "kdil" "d"]);

save('system_maturation_delay.mat', 'system')

load('system_maturation_delay.mat')

first_obs = 2;
dt = 5;
noise_level = .1;
seed = 1;

generator = Generator(system ...                                            % generator setup
                      , 'N', 10 ...                                        % number of cells
                      , 't', unique([0 5 10 20 dt:dt:200]) ...              % time grid
                      , 'error_std', noise_level ...                        % std of lognormal multiplicative errors
                      , 'D_mult', .25 ...                                   % covariance matrix s.t. D = diag(D_mult*beta)^2
                      , 'observed', first_obs:system.K ...                  % observed states labels (or indices)
                      );

rng(seed);
[data, ground_truth] = generator.generate();
plot(generator)
% data.beta = ground_truth.beta;


%% Estimate ----------------------------------------------------------------

methods = [];
methods = [methods "GMGTS"];
% methods = [methods "GTS"];

estimator = Estimator(system, data.traces, data.t, 2 ...                    % estimator setup
                      , 'Stages', 2 ...                                     % 0: smoothing only, 1: first stage only
                      , 'Methods', methods ...                              % GMGT, GTS, or both
                      , 'TestConvergence', true ...
                      );

estimator.estimate();

close all
plot(estimator, 'True', ground_truth ...
     , 'States', 1:6 ...
     , 'MaxCells', 7)

