%% Setup system ------------------------------------------------------------
                          
clearvars
close all

% model =  'model_repressilator_full.txt';
% system = ODEIQM(model, 'FixedParameters', ["DNAT" "kf" "Kd" "m1" "p1"]);
% 
% save('system_repressilator.mat', 'system')

load('system_repressilator.mat')

observed = [1:3 7:9];
seed = 1;

correlation = eye(6) + diag([.5 0 .5 0 .5], 1) + diag([.5 0 .5 0 .5], -1);
D = .01 * system.k0 .* system.k0' .* correlation;

generator = Generator(system ...                                            % generator setup
                      , 'N', 100 ...                                        % number of cells
                      , 't', [0:5:100] ...                                  % time grid
                      , 'error_std', .05 ...                               % std of lognormal multiplicative errors
                      , 'D', D ...                                   % variance scale
                      , 'observed', observed ...             % observed states labels (or indices)
                      ...%, 'lognormal', true ...
                      );

rng(seed);
[data, ground_truth] = generator.generate();
% plot(generator)


%% Estimate ----------------------------------------------------------------

methods = [];
methods = [methods "GMGTS"];
methods = [methods "GTS"];

estimator = Estimator(system, data ...                                      % estimator setup
                      , 'Stages', 2 ...                                     % 0: smoothing only, 1: first stage only
                      , 'Methods', methods ...                              % GMGT, GTS, or both
                      ...%, 'LogNormal', true ...
                      );

estimator.estimate();

close all
plot(estimator, 'True', ground_truth ...
     , 'States', 1:9 ...
     , 'MaxCells', 7)

