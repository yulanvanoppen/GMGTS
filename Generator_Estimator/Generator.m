classdef Generator < handle
    properties (SetAccess = private)
        system                                                              % nested object controlling ODE system
        settings = struct();                                                % data generation settings
        data = struct();                                                    % generated data
    end
    
    
    methods
        function obj = Generator(system, varargin)                      % Constructor
            default_N = 20;                                                 % number of cells
            default_t = [0 5 10 20:10:200];                                 % time grid
            default_error_std = .05;                                        % std of lognormal multiplicative errors
            default_init = system.x0' + 1e-4;                               % boundary conditions
            default_b = system.k0';                                         % initial estimate
            default_LogNormal = false;

            parser = inputParser();
            
            overlaps = @(x, S) ~isempty(intersect(x, S));
            addRequired(parser, 'system', @(x) isa(x, 'System') || isstring(string(x)) && numel(string(x)) == 1);
            addParameter(parser, 'N', default_N, @(x) isnumeric(x) && x >= 2);
            addParameter(parser, 't', default_t, @isnumeric);
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
                    default_D_mult = .01;
                    
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
            
            errors = normrnd(0, obj.settings.error_std, size(obj.data.original));
            obj.data.traces = obj.data.original .* (1 + errors);            % multiplicative measurement noise

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