%% Setup system ------------------------------------------------------------
                          
clearvars
close all

% name = 'bifunctional_TCS';
% model =  'model_bifunctional.txt';
% system = ODEIQM(name, model , 'FixedParameters', ["k1" "k2" "k4" "k6"]);
% 
% save('system_bifunctional_measurable.mat', 'system')

load('system_bifunctional_measurable.mat')

first_obs = 5;
dt = 10;
noise_level = .05;
seed = 1;

generator = Generator(system ...                                            % generator setup
                      , 'N', 100 ...                                        % number of cells
                      , 't', unique([0 2.5 5 dt:dt:100]) ...                % time grid
                      , 'error_std', noise_level ...                        % std of lognormal multiplicative errors
                      , 'D_mult', .25 ...                                   % variance scale
                      , 'observed', first_obs:system.K ...                                 % observed states labels (or indices)
                      , 'varying', 1:system.P ...                           % variable parameter labels (or indices)
                      );

rng(seed);
generator.generate();                                                       % generate data and package
% plot(generator);

generated = generator.data;
[data, ground_truth] = obfuscate(generated);
data.beta = ground_truth.beta;

[knots, penalized] = inflections(data)


%% Estimate ----------------------------------------------------------------

methods = [];
methods = [methods "GMGTS"];
% methods = [methods "GTS"];

% estimator = Estimator(data, system ...                                      % estimator setup
%                       , 'Stages', 2 ...                                     % 0: smoothing only, 1: first stage only
%                       , 'Methods', methods ...                              % GMGT, GTS, or both
%                       , 'Knots', [0 2.5 40] ...
%                       , 'PenalizedInterval', [5 100] ...
%                       );
                  
estimator = Estimator(data, system ...                                      % estimator setup
                      , 'Stages', 2 ...                                     % 0: smoothing only, 1: first stage only
                      , 'Methods', methods ...                              % GMGT, GTS, or both
                      , 'Knots', [0 2.5 linspace(20, 100, 4)] ...
                      , 'PenalizedInterval', penalized ...
                      , 'Lambda', 1e-12 ...
                      );

estimator.estimate();

% close all
plot(estimator, 'True', ground_truth ...
     , 'States', 1:6 ...
     , 'MaxCells', 7)

