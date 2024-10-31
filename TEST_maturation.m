%% Setup system ------------------------------------------------------------
                          
clearvars
close all

% name = 'maturation_fluorescence'; 
% model =  'model_maturation_onestep.txt';
% system = ODEIQM(name, model, 'FixedParameters', ["kr" "kdr" "kdil" "d"]);
% 
% save('system_maturation_delay.mat', 'system')

load('system_maturation_delay.mat')

first_obs = 1;
dt = 5;
noise_level = .1;
seed = 1;

generator = Generator(system ...                                            % generator setup
                      , 'N', 100 ...                                        % number of cells
                      , 't', unique([0 5 10 20 dt:dt:200]) ...                                  % time grid
                      , 'error_std', noise_level ...                               % std of lognormal multiplicative errors
                      , 'D_mult', .25 ...
                      , 'observed', first_obs:system.K ...                            % observed states labels (or indices)
                      , 'varying', 1:system.P ...                           % variable parameter labels (or indices)
                      );

rng(seed);
generator.generate();                                                       % generate data and package
% plot(generator);

generated = generator.data;
[data, ground_truth] = obfuscate(generated);
data.beta = ground_truth.beta;
data.original = ground_truth.original;
data.doriginal = ground_truth.doriginal;

% knots = placement(data)

%% Estimate ----------------------------------------------------------------

methods = [];
methods = [methods "GMGTS"];
% methods = [methods "GTS"];

estimator = Estimator(data, system ...                                      % estimator setup
                      , 'Stages', 2 ...                                     % 0: smoothing only, 1: first stage only
                      , 'Methods', methods ...                              % GMGT, GTS, or both
                      , 'InteractiveSmoothing', true ...
                      );

estimator.estimate();

% close all
plot(estimator, 'True', ground_truth ...
     , 'States', 1:6 ...
     , 'MaxCells', 7)
