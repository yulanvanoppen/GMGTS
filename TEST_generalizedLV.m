%% Setup system ------------------------------------------------------------
                          
clearvars
close all

% model =  'model_generalizedLV2.txt';
% system = ODEIQM(model, 'FixedParameters', strcat('r', string(1:16)));
% save('system_generalizedLV2.mat', 'system')

load('system_generalizedLV2.mat')

first_obs = 10;
noise_level = .05;
seed = 1;

generator = Generator(system ...                                            % generator setup
                      , 'N', 10 ...                                         % number of cells
                      , 't', 0:20 ...                                       % time grid
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

estimator = Estimator(system, data ...                                      % estimator setup
                      , 'Stages', 2 ...                                     % 0: smoothing only, 1: first stage only
                      , 'Methods', methods ...                              % GMGT, GTS, or both
                      , 'MaxIterationsFS', 50 ...
                      , 'ConvergenceTolFS', 1e-12 ...
                      , 'TestConvergence', true ...
                      );

estimator.estimate();

close all
plot(estimator, 'True', ground_truth ...
     , 'States', 9:16 ...
     , 'MaxCells', 7)

