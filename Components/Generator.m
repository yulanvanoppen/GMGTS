classdef Generator < handle
%GENERATOR handles data generation from ODE-based mixed-effect models
%
%   generator = GENERATOR(system, ...) instantiates a GENERATOR object to
%   draw parameters from a specified random effects distribution, integrate
%   ODE systems to obtain cell trajectories, and produce measurements by
%   perturbing the trajectories with measurement noise.
%
%   measurements = generator.generate() generates a struct containing
%   fields 'traces' (perturbed trajectories), 't' (measurement time points), 'observed' (observed state indices),
%   'init' (assumed initial conditions), 'T' (number of time points), 'L'
%   (number of observables), and 'N' (number of cells).
%
%   measurements = generator.generate(beta) uses the cell-specific 
%   (P-dimensional) parameter vectors in beta (NxP-dimensional) instead of
%   randomly drawn random effects (useful for boostrapping tests).
%
%   [measurements, ground_truth] = generator.generate(...) additionally
%   returns the complete ground_truth struct, which also contains the 
%   population parameters, the cell-specific parameters, and the unperturbed
%   cell trajectories.
%
%   plot(generator) plots generated trajectories and corresponding 
%   underlying gradients for each state.
%
%   CONSTRUCTOR NAME-VALUE ARGUMENTS
%     N - Number of cells to generate
%       20 (default) | positive integer
%
%     t - Measurement times
%       0:20:200 (default) | numeric vector
%
%     error_const - Additive noise standard deviation
%       0 (default) | positive scalar
%
%     error_std - Multiplicative noise standard deviation
%       .05 (default) | positive scalar
%
%     init - Initial conditions
%       numeric vector
%       (uses the initial conditions specified by system by default)
%
%     b - Random effects mean vector
%       numeric vector
%       (uses the nominal parameter vectors of system by default)
%
%     D_mult - Random effects (common) coefficient of variation
%       .1 (default) | positive scalar
%       (ignored when D is specified, see below)
%
%     D - Random effects covariance matrix
%       positive semidefinite matrix
%
%     observed - Observed state indices
%       1:system.K (default) | positive integer vector
%
%     lognormal - Use log-normal random effects distribution
%       false (default) | true
%       (if true, the log-mean Lb and log-covariance matrix LD are approximated
%       by moment matching: LD==log(1+D./(b'*b)) and Lb==log(b)-diag(LD)/2)
% 
%   See also GMGTS, SYSTEM, ESTIMATOR

    properties (SetAccess = private)
        system                                                              % nested object controlling ODE system
        settings = struct();                                                % data generation settings
        data = struct();                                                    % generated data
    end
    
    
    methods
        function obj = Generator(system, varargin)                      % Constructor
            default_N = 20;                                                 % number of cells
            default_t = 0:20:200;                                           % time grid
            default_error_const = 0;                                        % additive error standard deviation
            default_error_std = .05;                                        % std of 'multiplicative' errors
            default_init = system.x0' + 1e-8;                               % boundary conditions
            default_b = system.k0';                                         % initial estimate
            default_LogNormal = false;

            parser = inputParser();
            
            overlaps = @(x, S) ~isempty(intersect(x, S));
            addRequired(parser, 'system', @(x) isa(x, 'System') || isstring(string(x)) && numel(string(x)) == 1);
            addParameter(parser, 'N', default_N, @(x) isnumeric(x) && x >= 2);
            addParameter(parser, 't', default_t, @isnumeric);
            addParameter(parser, 'error_const', default_error_const, @(x) x > 0);
            addParameter(parser, 'error_std', default_error_std, @(x) x > 0);
            addParameter(parser, 'init', default_init, ...
                         @(x) isnumeric(x) && numel(x) == system.K);
            addParameter(parser, 'b', default_b, ...
                         @(x) isnumeric(x) && numel(x) == system.P);
            addParameter(parser, 'D_mult', [], @(x) isempty(x) || x > 0);
            addParameter(parser, 'D', [], ...
                         @(x) isempty(x) || issymmetric(x) && all(eig(x) > -1e-15) ...
                                                           && size(x, 1) == system.P);
            addParameter(parser, 'observed', [], ...
                         @(x) isempty(x) || ...
                              overlaps(string(x), string(1:system.K)) || ...
                              overlaps(string(x), string(system.states)) );
            addParameter(parser, 'lognormal', default_LogNormal, @islogical);
            
            parse(parser, system, varargin{:});
            
            obj.settings = parser.Results;
            
            if isstring(obj.settings.observed)
                observed = 1:system.K;
                obj.settings.observed = observed(ismember(system.states, ...
                                                          obj.settings.observed));
            end
            
            obj.system = system;                                            % ODE system handle
            if ~isa(system, 'System')                                       % process model file if provided
                system = ODEIQM(string(system), varargin{:});
            end
            obj.system = system;
            obj.add_system_defaults();                                      % complete with defaults
        end
        
        
        function add_system_defaults(obj)                               % Supply default settings
            default_D_mult = .1;                                            % variance multiplier scale
            default_observed = 1:obj.system.K;                              % observed state indices
            
            observed_labels = [];                                           % observed state labels
            
            switch obj.system.name
                case 'bifunctional_TCS'
                    default_D_mult = .25;
                    
                case 'maturation_fluorescence'
                    default_D_mult = .25;
                    
                case 'repressilator'
                    default_D_mult = .1;
                    
                case 'generalizedLV'
                    default_D_mult = .25;
            end
            
            if isempty(obj.settings.D_mult)
                D_mult = default_D_mult;
            else
                D_mult = obj.settings.D_mult;
            end
            default_D = diag((D_mult * obj.settings.b).^2);
            
            if ~isempty(observed_labels)
                default_observed = default_observed(ismember(obj.system.states, observed_labels));
            end
            
            if isempty(obj.settings.D), obj.settings.D = default_D; end
            if isempty(obj.settings.observed), obj.settings.observed = default_observed; end
        end
        
        
        function [measurements, ground_truth] = generate(obj, beta)     % Generate data
            obj.data.b = obj.settings.b;                                    % copy relevant settings
            obj.data.D = obj.settings.D;
            obj.data.init = obj.settings.init;
            obj.data.t = obj.settings.t;
            obj.data.observed = obj.settings.observed;
            obj.data.lognormal = obj.settings.lognormal;
            obj.data.T = length(obj.settings.t);
            obj.data.L = length(obj.settings.observed);
            obj.data.N = obj.settings.N;
            
            if nargin < 2                                                   % draw nonnegative random effects
                if obj.data.lognormal
                    obj.data.D = log(1 + obj.data.D ./ (obj.data.b'*obj.data.b));
                    obj.data.b = log(obj.data.b) - .5 * diag(obj.data.D)';
                    obj.data.beta = exp(mvnrnd(obj.data.b, obj.data.D, obj.data.N));
                else
                    obj.data.beta = abs(mvnrnd(obj.data.b, obj.data.D, obj.data.N));
                end
            else
                obj.data.beta = beta;
            end
                                  
            if obj.data.lognormal                                          % numerically integrate ODE system
                obj.data.moriginal = obj.system.integrate(exp(obj.data.b), obj.settings);
                obj.data.mdoriginal = obj.system.rhs(obj.data.moriginal, obj.data.t, exp(obj.data.b));
            else
                obj.data.moriginal = obj.system.integrate(obj.data.b, obj.settings);
                obj.data.mdoriginal = obj.system.rhs(obj.data.moriginal, obj.data.t, obj.data.b);
            end
                                                                            % corresponding cell-specific traces
            obj.data.original = obj.system.integrate(obj.data.beta, obj.settings);
            obj.data.doriginal = obj.system.rhs(obj.data.original, obj.data.t, obj.data.beta);
            
            add = normrnd(0, obj.settings.error_const, size(obj.data.original));
            mul = normrnd(0, obj.settings.error_std, size(obj.data.original));
            obj.data.traces = add + obj.data.original .* (1+mul);           % additive/multiplicative measurement noise
            obj.data.noisevar = obj.settings.error_const^2 + obj.settings.error_std^2 * obj.data.original.^2;

            [measurements, ground_truth] = obj.separate();
        end


        function [measurements, ground_truth] = separate(obj)           % Separate simulated measurements from ground truth
            ground_truth = obj.data;
            measurements = struct('traces', obj.data.traces(:, obj.data.observed, :), ...
                                  't', obj.data.t, 'observed', obj.data.observed, ...
                                  'init', obj.data.init, 'T', obj.data.T, 'L', obj.data.L, 'N', obj.data.N);
        end
    
        
        
        function plot(obj)
            if ~isfield(obj.data, 'traces')
                obj.generate();
            end
            
            set(groot, 'defaultAxesTickLabelInterpreter','latex');
            set(groot, 'defaultLegendInterpreter','latex');
            set(groot,'defaultTextInterpreter','latex');
            
            screen = get(0, 'ScreenSize');
            figure('position', [10, 70, min(1740, screen(3)-180), min(900, screen(4)-180)])
            tl = tiledlayout(ceil(obj.system.K / 5), min(obj.system.K, 5));
            n_cells = min(obj.data.N, 20);
            title(tl, sprintf('Generated data (%d/%d shown)', n_cells, obj.data.N))
            
            for k = 1:obj.system.K
                nexttile(k)
                title(obj.system.states(k))
                hold on
                yyaxis left
                ylabel('measurement')
                h1 = plot(obj.data.t, reshape(obj.data.traces(:, k, 1:n_cells), ...
                                              obj.data.T, n_cells), '-');
                h2 = plot(obj.data.t, reshape(obj.data.original(:, k, 1:n_cells), ...
                                              obj.data.T, n_cells), '--');
                yyaxis right
                ylabel('gradient')
                h3 = plot(obj.data.t, reshape(obj.data.doriginal(:, k, 1:n_cells), ...
                                              obj.data.T, n_cells), '--');
                xlabel('time')
                
                set(h1, {'color'}, num2cell(winter(n_cells), 2));
                set(h2, {'color'}, num2cell(winter(n_cells), 2));
                set(h3, {'color'}, num2cell(autumn(n_cells), 2));
            end
        end
    end
end