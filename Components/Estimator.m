classdef Estimator < handle
%ESTIMATOR manages the stages of the GTS and GMGTS methods.
%
%   CONSTRUCTOR
%   estimator = ESTIMATOR(system, data, ...) instantiates an ESTIMATOR
%   object to infer random effect distributions of the specified
%   K-dimensional (ODE) system (see the System class) from measurements
%   given in data. Here data is either a TxLxN-dimensional array with
%   measurements at T time points for L observables and N cells
%   (individuals), or a 1x1 struct with its measurements stored in a field
%   named 'y' or 'traces' . This struct may additionally contain fields
%   't', a (T-dimensional) array of measurement time points, 'observed', an
%   (L-dimensional) vector of indices of observed states (with respect to
%   the specified system), and 'init', a (K-dimensional) vector of initial
%   conditions. Default values are imputed for any omitted additional
%   fields.
%
%   estimator = ESTIMATOR(system, traces, t, observed, init, ...) provides
%   an alternative constructor with (optional) positional arguments.
%   Omitting t, observed, or init causes default values to be used for the
%   omitted and all subsequent arguments. The defaults are:
%       t = 0:size(data, 1)-1
%       observed = 1:size(data, 2)
%       init = 1e-8*ones(1, system.K)
%
%   OBJECT METHODS
%       varargout = estimator.estimate() carries out estimation for the
%       specified stages and method(s), returning two arguments when both
%       both methods are used.
%
%       varargout = estimator.estimate(silent) suppresses console prompts 
%       during estimation if silent==true.
%       
%       plot(estimator, ...) produces relevant plots for the executed stages
%       and methods.
%
%       plot(estimator, True=ground_truth, ...) additionally plots the
%       data-generating parameters and trajectories contained in the
%       ground_truth struct for reference; see the Generator class.
%
%       plot(estimator, States=states, ...) restricts the plots to the
%       specified states, which may consist of indices or state labels.
%
%       plot(estimator, Parameters=parameters, ...) restricts the plots to
%       the specified parameters, which may consist of indices or labels.

%       plot(estimator, MaxCells=N, ...) only plots data and predictions for
%       the first N cells.
% 
%   CONSTRUCTOR NAME-VALUE ARGUMENTS
%     Stages - GMGTS and GTS stages to execure
%       2 (default) | 0 | 1
%       (GMGTS either executes smoothing only (0), smoothing
%       and first-stage estimates (1), or everything including the
%       second-stage estimates (2); GTS executes the first-stage only (0,
%       1) or everything (2))
%
%     Methods - Frameworks to use
%       "GMGTS" (default) | ["GMGTS" "GTS"] | "GTS"
%        
%     AutoKnots - Use automatic B-spline knot placement heuristic
%       true (default) | false
%
%     Knots - Knot locations used for B-spline smoothing
%       numeric vector | cell array
%       (either a single numeric vector of knot locations for all observed
%       states, or a cell array of state-specific knot location vectors)
%
%     InteractiveSmoothing - Use the interactive app to smooth measurements
%       true (default) | false
%
%     LB - Parameter space lower bounds 
%       numeric vector
%       (.25 times the nominal parameter values specified by system by default)
%
%     UB - Parameter space upper bounds
%       numeric vector
%       (4 times the nominal parameter value specified by system by default)
%
%     PositiveState - Force smoothed/predicted state positivity
%       true (default) | false
%
%     TimePoints - Time point grid for first-stage optimization
%       numeric vector
%       (10 equidistant intervals in the measurement interval by default)
%
%     MaxIterationsSM - Maximum number of smoothing iterations
%       20 (default) | positive integer
%
%     ConvergenceTolSM - Convergence step tolerance for smoothing iterations
%       1e-3 (default) | positive scalar
%
%     NMultiStartFS - Number of multistarts for first-stage initialization
%       10 (default) | positive integer
%
%     MaxIterationsFS - Maximum number of first-stage iterations
%       5 (default) | positive integer
%
%     ConvergenceTolFS - Convergence step tolerance for first-stage iterations
%       2e-3 (default) | positive scalar
%
%     MaxIterationsSS - Maximum number of second-stage iterations
%       10 (default) | positive integer
%
%     ConvergenceTolSS - Convergence step tolerance for second-stage iterations
%       1e-3 (default) | positive scalar
%
%     LogNormal - Infer log-normal random effect distribution
%       false (default) | true
%
%     Prior - Include parameter prior (same for each cell)
%       struct('mean', 0, 'prec', 0) (default) | struct
%       (the struct should have fields 'mean' specifying the prior mean,
%       and either 'prec', the prior precision matrix, 'sd', the prior
%       component-wise standard deviations, or 'cv', the prior component-
%       wise coefficients of variation; if multiple fields indicating
%       variation are included, only one is considered with precedence 'cv',
%       'sd', 'prec' from high to low)
%
%   See also GMGTS, SYSTEM, GENERATOR
    
    properties (SetAccess = private)
        data                                                                % measurements struct
        system                                                              % ODE system object
        
        stages                                                              % stages to execute
        method                                                              % method(s) to conduct inference
        autoknots                                                           % logical for placement heuristoic
        knots                                                               % B-spline knots for each state
        interactive                                                         % logical for interactive smoothing app use
        
        lb                                                                  % parameter space lower bounds
        ub                                                                  % parameter space upper bounds
        positive                                                            % logical for forced state positivity
        t_fs                                                                % first stage optimization grid (GMGTS)
        
        nmultistart                                                         % #starting points for first stage initialization
        niterSM                                                             % #iterations for smoothing
        tolSM                                                               % convergence tolerance for smoothing
        niterFS                                                             % #iterations for the first stage
        tolFS                                                               % convergence tolerance for the first stage
        niterSS                                                             % #iterations for the second stage 
        tolSS                                                               % convergence tolerance for the second stage
        
        testconv                                                            % pipeline to approximate FS basins of attraction
        lognormal                                                           % whether to infer lognormal RE distribution
        prior                                                               % parameter prior distribution
        
        GMGTS_settings                                                      % GMGTS settings struct
        GMGTS_smoother                                                      % GMGTS smoothing object
        GMGTS_first_stage                                                   % GMGTS first stage optimization object
        GMGTS_second_stage                                                  % GMGTS second stage inference object
        GMGTS_conv_test                                                     % GMGTS basin of attraction approximator
        results_GMGTS                                                       % collect GMGTS results
        
        GTS_settings                                                        % GTS settings struct
        GTS_first_stage                                                     % GTS first stage optimization object
        GTS_second_stage                                                    % GTS second stage inference object
        results_GTS                                                         % collect GTS results
    end
    
    
    methods
        %% Constructor -----------------------------------------------------
        function obj = Estimator(system, data, varargin)
            [system, data] = obj.parse_initial(system, data, varargin{:});
            
            default_t = 0:size(data, 1)-1;
            default_observed = 1:size(data, 2);
            default_init = 1e-8 * ones(1, system.K);
            
            default_Stages = 2;
            default_Methods = "GMGTS";
            
            initial = system.k0';
            default_AutoKnots = true;
            default_Knots = repmat({linspace(data.t(1), data.t(end), round((data.T-1)/2)+1)}, ...
                                   1, length(data.observed));
            default_InteractiveSmoothing = false;
            
            default_LB = .25 * initial;
            default_UB = 4 .* initial + .0001 * mean(initial);
            default_PositiveStates = true;
            default_TimePoints = data.t(1) + (0:.1:1) * range(data.t);
            
            default_MaxIterationsSM = 20;
            default_ConvergenceTolSM = 1e-3;
            default_NMultiStartFS = 10;
            default_MaxIterationsFS = 5;
            default_ConvergenceTolFS = 2e-3;
            default_MaxIterationsSS = 10;
            default_ConvergenceTolSS = 1e-3;
            
            default_TestConvergence = false;
            default_LogNormal = false;
            default_Prior = struct('mean', 0, 'prec', 0);
            
            parser = inputParser;
            parser.KeepUnmatched = true;
            addRequired(parser, 'system', @(x) isa(x, 'System') || isstring(string(x)) && numel(string(x)) == 1);
            addRequired(parser, 'data', @(x) isstruct(x) || isnumeric(x) && ndims(x) == 3 && size(x, 1) > 1);
            addOptional(parser, 't', default_t, @(x) isnumeric(x) && numel(unique(x)) == data.T);
            addOptional(parser, 'observed', default_observed, @(x) isnumeric(x) && numel(unique(x)) == data.L);
            addOptional(parser, 'init', default_init, @(x) isnumeric(x) && numel(x) == system.K);
            
            addParameter(parser, 'Stages', default_Stages, @(x) ismember(x, 0:2));
            addParameter(parser, 'Methods', default_Methods, ...
                         @(x) all(ismember(upper(string(x)), ["GMGTS" "GTS"])));
            addParameter(parser, 'AutoKnots', default_AutoKnots, @islogical);
            addParameter(parser, 'Knots', default_Knots, @(x) (iscell(x) && length(x) == data.L ...
                                                               && all(cellfun(@isnumeric, x))) ...
                                                           || isnumeric(x));
            addParameter(parser, 'InteractiveSmoothing', default_InteractiveSmoothing, @islogical);

            addParameter(parser, 'LB', default_LB, @(x) all(x < initial));
            addParameter(parser, 'UB', default_UB, @(x) all(x > initial));
            addParameter(parser, 'PositiveStates', default_PositiveStates, @islogical);
            addParameter(parser, 'TimePoints', default_TimePoints, @(x) all(data.t(1) <= x & x <= data.t(end)));
            
            addParameter(parser, 'MaxIterationsSM', default_MaxIterationsSM, @isscalar);
            addParameter(parser, 'ConvergenceTolSM', default_ConvergenceTolSM, @isscalar);
            addParameter(parser, 'NMultiStartFS', default_NMultiStartFS, @isscalar);
            addParameter(parser, 'MaxIterationsFS', default_MaxIterationsFS, @isscalar);
            addParameter(parser, 'ConvergenceTolFS', default_ConvergenceTolFS, @isscalar);
            addParameter(parser, 'MaxIterationsSS', default_MaxIterationsSS, @isscalar);
            addParameter(parser, 'ConvergenceTolSS', default_ConvergenceTolSS, @isscalar);
            
            addParameter(parser, 'TestConvergence', default_TestConvergence, @islogical);
            addParameter(parser, 'LogNormal', default_LogNormal, @islogical)
            addParameter(parser, 'Prior', default_Prior, @(x) isfield(x, 'mean') && numel(x.mean) == system.P ...
                                                           && (isfield(x, 'prec') && all(size(x.prec) == system.P) ...
                                                               && issymmetric(x.prec) && all(eig(x.prec) >= 0) ...
                                                            || isfield(x, 'sd') && numel(x.sd) == system.P ...
                                                               && all(x.sd > 0) ...
                                                            || isfield(x, 'cv') && numel(x.cv) == system.P ...
                                                               && all(x.cv > 0)));
            parse(parser, system, data, varargin{:});
            [obj.system, obj.data] = obj.parse_initial(parser.Results.system, parser.Results.data, varargin{:});
            obj.parse_parameters(parser);
            
            obj.data.T_fine = 81;
            obj.data.t_fine = linspace(obj.data.t(1), obj.data.t(end), obj.data.T_fine);
            
            if ismember("GMGTS", obj.method), obj.constructor_GMGTS; end
            if ismember("GTS", obj.method), obj.constructor_GTS; end
        end

                                                                        % Basic parsing to allow conditional input validations
        function [system, data] = parse_initial(~, system, data, varargin)  
            if ~isa(system, 'System')                                       % process model file if provided
                namevalue = arrayfun(@(idx) iscellstr(varargin(idx)) || isstring(varargin{idx}), 1:length(varargin));
                first_namevalue = find(namevalue, 1);
                if isempty(first_namevalue), first_namevalue = length(varargin)+1; end
                system = System(string(system), varargin{first_namevalue:end});
            end
            if isstruct(data)                                               % default any missing fields
                if ~isfield(data, 'traces') && isfield(data, 'y'), data.traces = data.y; end
                if ~isfield(data, 't'), data.t = 0:size(data.traces, 1)-1; end
                if ~isfield(data, 'observed'), data.observed = 1:size(data.traces, 2); end
                if ~isfield(data, 'init'), data.init = system.x0' + 1e-4; end
            else                                                            % components provided separately
                traces = data;                                              % array with measurements instead of struct
                if ~iscellstr(varargin(1)) && ~isstring(varargin{1})        % recursively check if optional or Name/Value
                    t = sort(unique(reshape(varargin{1}, 1, [])));
                    if ~iscellstr(varargin(2)) && ~isstring(varargin{2})
                        observed = sort(unique(reshape(varargin{2}, 1, [])));
                        if ~iscellstr(varargin(3)) && ~isstring(varargin{3})
                            init = reshape(varargin{3}, 1, []);
                        else                                                % substitute defaults instead
                            init = 1e-4 * ones(1, system.K);
                        end
                    else
                        observed = 1:size(traces, 2);
                        init = 1e-4 * ones(1, system.K);
                    end
                else
                    t = 0:size(traces, 1)-1;
                    observed = 1:size(traces, 2);
                    init = 1e-4 * ones(1, system.K);
                end                                                         % compile into struct
                data = struct('traces', traces, 't', t, 'observed', observed, 'init', init);
            end
            [data.T, data.L, data.N] = size(data.traces);                   % include dimensions for notational convenience
        end


        function parse_parameters(obj, parser)                          % Parse Name/Value constructor arguments
            obj.stages = parser.Results.Stages;
            obj.method = string(parser.Results.Methods);
            obj.autoknots = parser.Results.AutoKnots;
            obj.knots = cell(1, obj.data.L);
            parsed_knots = parser.Results.Knots;
            if ~iscell(parsed_knots)
                parsed_knots = repmat({parsed_knots}, 1, obj.data.L);
            end
            for state = 1:length(obj.data.observed)
                truncated = max(obj.data.t(1), min(obj.data.t(end), parsed_knots{state}));
                arranged = sort(unique([obj.data.t(1) reshape(truncated, 1, []) obj.data.t(end)]));
                obj.knots{state} = (arranged - obj.data.t(1)) / range(obj.data.t);
            end
            obj.autoknots = parser.Results.AutoKnots && ismember("Knots", string(parser.UsingDefaults));
            obj.interactive = parser.Results.InteractiveSmoothing;
            
            obj.lb = parser.Results.LB;
            obj.ub = parser.Results.UB;
            obj.positive = parser.Results.PositiveStates;
            obj.t_fs = sort(unique(parser.Results.TimePoints));
            obj.niterSM = max(1, round(parser.Results.MaxIterationsSM));
            obj.tolSM = max(1e-12, parser.Results.ConvergenceTolSM);
            obj.nmultistart = max(1, round(parser.Results.NMultiStartFS));
            obj.niterFS = max(1, round(parser.Results.MaxIterationsFS));
            obj.tolFS = max(1e-12, parser.Results.ConvergenceTolFS);
            obj.niterSS = max(1, round(parser.Results.MaxIterationsSS));
            obj.tolSS = max(1e-12, parser.Results.ConvergenceTolSS);
            
            obj.testconv = parser.Results.TestConvergence;
            obj.lognormal = parser.Results.LogNormal;
            obj.prior = parser.Results.Prior;
            obj.prior.mean = flatten(obj.prior.mean);
            if isfield(obj.prior, 'cv')
                obj.prior.cv = flatten(obj.prior.cv);
                obj.prior.mean(obj.prior.mean == 0 & isinf(obj.prior.cv)) = 1;
                obj.prior.sd = obj.prior.cv .* obj.prior.mean;
            end
            if isfield(obj.prior, 'sd')
                obj.prior.sd = flatten(obj.prior.sd);
                warning('off','MATLAB:singularMatrix')
                obj.prior.prec = inv(diag(obj.prior.sd.^2));
                warning('on','MATLAB:singularMatrix')
            end
            if ~isfield(obj.prior, 'mult'), obj.prior.mult = 1; end
        end
        
        
        function constructor_GMGTS(obj)
            weights = ones(length(obj.t_fs), obj.system.K);                 % omit interval ends
            weights([1 end], :) = 0;
            
            obj.GMGTS_settings.sm = struct('order', 4, 'autoknots', obj.autoknots, 'knots', {obj.knots}, ...
                                           'positive', obj.positive, 't_fs', obj.t_fs, 'niter', obj.niterSM, ...
                                           'tol', obj.tolSM, 'interactive', obj.interactive);
            obj.GMGTS_settings.fs = struct('lb', obj.lb, 'ub', obj.ub, 'positive', obj.positive, 't_fs', obj.t_fs, ...
                                           'weights', weights, 'niter', obj.niterFS, 'tol', obj.tolFS, ...
                                           'nstart', obj.nmultistart, 'lognormal', obj.lognormal, 'prior', obj.prior);
            obj.GMGTS_settings.ss = struct('positive', obj.positive, 'niter', obj.niterSS, 'tol', obj.tolSS, ...
                                           'lognormal', obj.lognormal);
        end
        
        
        function constructor_GTS(obj)
            optindices = 1:length(obj.data.t);                              % time point indices to restrict opt to
            
            obj.GTS_settings.fs = struct('lb', obj.lb, 'ub', obj.ub, 'optindices', optindices, ...
                                         'positive', obj.positive, 'niter', obj.niterFS, 'tol', obj.tolFS, ...
                                         'nstart', 10*obj.nmultistart, 'prior', obj.prior, ...
                                         'lognormal', obj.lognormal);
            obj.GTS_settings.ss = struct('positive', obj.positive, 'niter', obj.niterSS, 'tol', obj.tolSS, ...
                                         'lognormal', obj.lognormal);
        end
        
        
        %% Estimation ------------------------------------------------------
        function varargout = estimate(obj, silent)
            if nargin < 2, silent = false; end
            varargout = cell(1, 0);
            idx = 1;
            
            if ismember("GMGTS", obj.method)
                obj.estimate_GMGTS(silent);
                
                out = struct('t', obj.results_GMGTS.t_fine, 'smoothed', obj.results_GMGTS.smoothed_fine, ...
                             'dsmoothed', obj.results_GMGTS.dsmoothed_fine);
                if obj.stages >= 1
                    out.fitted = obj.results_GMGTS.fitted_fs_fine;
                    out.dfitted = obj.results_GMGTS.fitted_fs_fine;
                    out.beta = obj.results_GMGTS.beta_fs;
                    out.uncertainties = obj.results_GMGTS.varbeta;
                end
                if obj.stages == 2
                    out.b = obj.results_GMGTS.b_est;
                    out.D = obj.results_GMGTS.D_est;
                end
                
                varargout{idx} = out;
                idx = idx + 1;
            end
            
            if ismember("GTS", obj.method)
                obj.estimate_GTS(silent);
                
                out = struct('t', obj.results_GTS.t_fine, 'smoothed', [], 'dsmoothed', [], ...
                             'fitted', obj.results_GTS.fitted_fs_fine, ...
                             'dfitted', obj.results_GTS.dfitted_fs_fine, 'beta', obj.results_GTS.beta_fs, ...
                             'uncertainties',  obj.results_GTS.varbeta);
                if obj.stages == 2
                    out.b = obj.results_GTS.b_est;
                    out.D = obj.results_GTS.D_est;
                end
                    
                varargout{idx} = out;
            end
        end
        
        
        function estimate_GMGTS(obj, silent)
            ws = warning('error', 'MATLAB:nearlySingularMatrix');
            
            obj.GMGTS_smoother = Smoother(obj.data, obj.system, obj.GMGTS_settings.sm);
            obj.results_GMGTS = obj.GMGTS_smoother.smooth();
            toc_sm_GMGTS = obj.results_GMGTS.toc_sm;
            if ~silent, toc_sm_GMGTS, end
            toc_fs_GMGTS = 0; toc_ss_GMGTS = 0;
            
            if obj.stages >= 1
                tic
                obj.GMGTS_first_stage = FirstStage(obj.results_GMGTS, obj.system, ...
                                                   obj.GMGTS_settings.fs);
                obj.results_GMGTS = obj.GMGTS_first_stage.optimize();
                if ~silent, toc_fs_GMGTS = toc, else, toc_fs_GMGTS = toc; end
            end
            if obj.stages == 2
                tic
                obj.GMGTS_second_stage = SecondStage(obj.results_GMGTS, obj.system, obj.GMGTS_settings.ss);
                obj.results_GMGTS = obj.GMGTS_second_stage.optimize();
                if ~silent, toc_ss_GMGTS = toc, else, toc_ss_GMGTS = toc; end
            end
            
            if ~silent
                fprintf('total time (GMGTS): %d seconds\n', toc_sm_GMGTS + toc_fs_GMGTS + toc_ss_GMGTS)
            end
            
            obj.results_GMGTS.time = [toc_sm_GMGTS toc_fs_GMGTS toc_ss_GMGTS];
            
            if obj.testconv
                obj.GMGTS_conv_test = ConvTest(obj.results_GMGTS, obj.system, obj.GMGTS_settings.fs);
                obj.results_GMGTS = obj.GMGTS_conv_test.estimate();
            end
            
            warning(ws);
        end
        
        
        function estimate_GTS(obj, silent)
            tic
            obj.GTS_first_stage = FirstStageGTS(obj.data, obj.system, obj.GTS_settings.fs);
            obj.results_GTS = obj.GTS_first_stage.optimize();
            if ~silent, toc_fs_GTS = toc, else, toc_fs_GTS = toc; end
            toc_ss_GTS = 0;

            if obj.stages == 2
                tic
                obj.GTS_second_stage = SecondStage(obj.results_GTS, obj.system, ...
                                                   obj.GTS_settings.ss);
                obj.results_GTS = obj.GTS_second_stage.optimize();
                if ~silent, toc_ss_GTS = toc, else, toc_ss_GTS = toc; end
            end
            
            if ~silent
                fprintf('total time (GTS): %d seconds\n', toc_fs_GTS + toc_ss_GTS)
            end
            
            obj.results_GTS.time = [toc_fs_GTS toc_ss_GTS];
        end
        
        
        function L = loglik(obj, nrep, ~)
            if nargin < 3, selected_method = "GMGTS"; else, selected_method = "GTS"; end
            if ~isfield(obj.results_GMGTS, 'b_est'), L = -Inf; return, end
            
            L = zeros(1, nrep);
            if selected_method == "GMGTS"
                b = obj.results_GMGTS.b_est;
                D = obj.results_GMGTS.D_est;
                var = @(trace) obj.GMGTS_smoother.theta_fs(1) + obj.GMGTS_smoother.theta_fs(2) * trace.^2;
            else
                b = obj.results_GTS.b_est;
                D = obj.results_GTS.D_est;
                var = @(trace) obj.GTS_first_stage.theta_fs(1) + obj.GTS_first_stage.theta_fs(2) * trace.^2;
            end
            
            logfs = @(trace, i) flatten(sum(-1/2 * (obj.data.traces(:, :, i) - trace).^2 ./ var(trace) - 1/2 * log(var(trace)), [1 2]));
            
            for rep = 1:nrep
                fprintf('%d ', rep);
                M = 1000;
                Lis = zeros(1, obj.data.N);

                for i = 1:obj.data.N
%                     fprintf('%d ', i);
                    beta_i = obj.results_GMGTS.beta_fs(i, :, :);

                    sample_parameters = mvnrnd(beta_i, D/9, M);
                    sample_traces = obj.system.integrate(sample_parameters, obj.data, obj.data.t, 1e-2);
                    logps = log(mvnpdf(sample_parameters, b, D));
                    
                    div = 1;
                    while isinf(quantile(logps, .9))
                        div = div + 1;
                        sample_parameters = mvnrnd((beta_i+(div-1)*b)/div, div*D/9, M);
                        sample_traces = obj.system.integrate(sample_parameters, obj.data, obj.data.t, 1e-2);
                        logps = log(mvnpdf(sample_parameters, b, D));
                    end
                    
                    logqs = log(mvnpdf(sample_parameters, (beta_i+(div-1)*b)/div, div*D/9));

                    terms = logfs(sample_traces(:, obj.data.observed, :), i) + logps - logqs;
                    reg = -max(terms);
                    Lis(i) = log(sum(exp(terms + reg), 'all')/M) - reg;
                end
%                 fprintf('\n %4f\n', sum(Lis));
                L(rep) = sum(Lis);
            end
            fprintf('\n Mean: %.2f SD: %.2f\n', mean(L), std(L));
        end
        
        
        %% Plotting --------------------------------------------------------
        function plot(obj, varargin)
            set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
            set(groot, 'defaultLegendInterpreter', 'latex');
            set(groot,'defaultTextInterpreter', 'latex');
            
            parser = inputParser;
            parser.KeepUnmatched = true;
            default_States = 1:obj.system.K;
            default_Parameters = 1:obj.system.P;
            default_MaxCells = 5;
            
            overlaps = @(x, S) ~isempty(intersect(string(x), string(S)));
            
            addParameter(parser, 'True', struct, @(x) isstruct(x));
            addParameter(parser, 'States', default_States, ...
                         @(x) overlaps(x, obj.data.observed) || ...
                              overlaps(x, obj.system.states));
            addParameter(parser, 'Parameters', default_Parameters, ...
                         @(x) overlaps(x, 1:obj.system.P) || ...
                              overlaps(x, obj.system.parameters_variable));
            addParameter(parser, 'MaxCells', default_MaxCells, @(x) x >= 1);
            parse(parser, varargin{:});
            plot_settings = rmfield(parser.Results, 'True');
            
            if isstring(plot_settings.States)
                all_states = 1:obj.system.K;
                plot_settings.States = all_states(ismember(obj.system.states, ...
                                                           plot_settings.States));
            else
                plot_settings.States = intersect(plot_settings.States, 1:obj.system.K);
            end
            if isstring(plot_settings.Parameters)
                all_parameters = 1:obj.system.P;
                plot_settings.Parameters = all_parameters(ismember(obj.system.parameters_variable, ...
                                                                   plot_settings.Parameters));
            else
                plot_settings.Parameters = intersect(plot_settings.Parameters, 1:obj.system.P);
            end
            
            if any(obj.method == "GMGTS")
                if ~isfield(obj.results_GMGTS, 'smoothed'), obj.estimate(); end
                obj.plot_smoothing(parser.Results.True, plot_settings)
            end
            
            if obj.stages >= 1
                obj.plot_individual_fits(parser.Results.True, plot_settings)
            end
            
            if obj.stages == 2
                obj.plot_population_fit(parser.Results.True, plot_settings)
                obj.plot_mixed_effects(parser.Results.True, plot_settings)
            end
        end
        
        
        function plot_smoothing(obj, truth, plot_settings)
            [states, ~, ind] = intersect(sort(plot_settings.States), obj.data.observed);
            n_cells = min(obj.data.N, plot_settings.MaxCells);
            
            screen = get(0, 'ScreenSize');
            figure('position', [25, 55, min(1740, screen(3)-180), min(900, screen(4)-180)])
            tl = tiledlayout(2, length(states));
            title(tl, sprintf('Smoothing (%d/%d shown)', n_cells, obj.data.N))
            
            for k = 1:length(states)
                nexttile(k)
                state_idx = obj.data.observed(ind(k));
                title(obj.system.states{state_idx})
                hold on
                ylabel('measurement')
                h1 = plot(obj.data.t, reshape(obj.data.traces(:, ind(k), 1:n_cells), ...
                          obj.data.T, n_cells), '-');
%                 h2 = plot(obj.data.t, reshape(obj.results_GMGTS.smoothed(:, ind(k), 1:n_cells), ...
%                           obj.data.T, n_cells), '--');
                h2 = plot(linspace(obj.data.t(1), obj.data.t(end), 81), ...
                          reshape(obj.results_GMGTS.smoothed_fine(:, ind(k), 1:n_cells), ...
                          [], n_cells), '--');
                if k == 1, legend([h1(1) h2(1)], 'data', 'smoothing'), end
                if ~isempty(fieldnames(truth))
                    h3 = plot(obj.data.t, reshape(truth.original(:, state_idx, 1:n_cells), ...
                              obj.data.T, n_cells), 'o');
                    if k == 1, legend([h1(1) h2(1) h3(1)], 'data', 'smoothing', 'truth'), end
                end
                xlabel('time')
                
                nexttile(length(states) + k)
                title(['$d$ \hspace*{-.3em}' obj.system.states{state_idx} '\hspace*{-.3em} $/dt$'])
                hold on
                ylabel('gradient')
%                 h4 = plot(obj.data.t, reshape(obj.results_GMGTS.dsmoothed(:, ind(k), 1:n_cells), ...
%                           obj.data.T, n_cells), '--');
                h4 = plot(linspace(obj.data.t(1), obj.data.t(end), 81), ...
                          reshape(obj.results_GMGTS.dsmoothed_fine(:, ind(k), 1:n_cells), ...
                          [], n_cells), '--');
                if k == 1, legend(h4(1), 'smoothing'), end
                if ~isempty(fieldnames(truth))
                    h5 = plot(obj.data.t, reshape(truth.doriginal(:, state_idx, 1:n_cells), ...
                              obj.data.T, n_cells), 'o');
                    if k == 1, legend([h4(1) h5(1)], 'smoothing', 'truth'), end
                end
                xlabel('time')
                
                set(h1, {'color'}, num2cell(parula(n_cells), 2));
                set(h2, {'color'}, num2cell(parula(n_cells), 2));
                set(h4, {'color'}, num2cell(parula(n_cells), 2));
                if ~isempty(fieldnames(truth))
                    set(h3, {'color'}, num2cell(parula(n_cells), 2));
                    set(h5, {'color'}, num2cell(parula(n_cells), 2));
                end
            end
        end
        
        
        function plot_individual_fits(obj, truth, plot_settings)
            states = sort(plot_settings.States);
            n_cells = min(obj.data.N, plot_settings.MaxCells);
            
            screen = get(0, 'ScreenSize');
            figure('position', [40, 40, min(1740, screen(3)-180), min(900, screen(4)-180)])
            tl = tiledlayout(1 + any(obj.method == "GMGTS"), length(states));
            title(tl, sprintf('Fitted individual trajectories (%d/%d shown)', n_cells, obj.data.N))
            
            col_scheme = num2cell(parula(n_cells), 2);
            for k = 1:length(states)
                nexttile(k)
                title(obj.system.states{states(k)})
                hold on
                if ismember(states(k), obj.data.observed)
                    obs_ind = find(states(k)==obj.data.observed);
                    h = plot(obj.data.t, reshape(obj.data.traces(:, obs_ind, 1:n_cells), ...
                                                  obj.data.T, n_cells), '-');
                    set(h, {'color'}, col_scheme);
                end
                xlabel('time')
                ylabel('measurement')
                if ismember(states(k), obj.data.observed) && states(k) == min(intersect(states, obj.data.observed))
                    legends = legend(h(1), 'data', 'AutoUpdate', 'off');
                end

                if any(obj.method == "GMGTS")
                    nexttile(length(states) + k)
                    title(['$d$ \hspace*{-.3em}' obj.system.states{k} '\hspace*{-.3em} $/dt$'])
                    hold on
                    if ismember(states(k), obj.data.observed)
                        h = plot(linspace(obj.data.t(1), obj.data.t(end), 81), ...
                                 reshape(obj.results_GMGTS.dsmoothed_fine(:, obs_ind, 1:n_cells), ...
                                         [], n_cells), '--');
                        set(h, {'color'}, col_scheme);
                    end
                    xlabel('time')
                    ylabel('gradient')
                    if ismember(states(k), obj.data.observed) && states(k) == min(intersect(states, obj.data.observed))
                        legends = [legends legend(h(1), 'smoothing', 'AutoUpdate', 'off')]; %#ok<AGROW>
                    end
                end
            end
            
            if numel(obj.method) == 1
                if obj.method == "GMGTS"
                    obj.add_individual_plots(obj.results_GMGTS.t_fine, obj.results_GMGTS.fitted_fs_fine, 'predictions', ...
                                             plot_settings, legends, ':', ...
                                             obj.results_GMGTS.dfitted_fs_fine, 'predictions')
                else
                    obj.add_individual_plots(obj.results_GTS.t, obj.results_GTS.fitted_fs, 'predictions', ...
                                             plot_settings, legends, '-.')
                end
            else
                obj.add_individual_plots(obj.results_GMGTS.t_fine, obj.results_GMGTS.fitted_fs_fine, 'predictions (GMGTS)', ...
                                         plot_settings, legends, ':', ...
                                         obj.results_GMGTS.dfitted_fs_fine, 'predictions (GMGTS)')
                obj.add_individual_plots(obj.results_GTS.t, obj.results_GTS.fitted_fs, 'predictions (GTS)', ...
                                         plot_settings, legends, '-.')
            end
            if ~isempty(fieldnames(truth))
                obj.add_individual_plots(obj.data.t, truth.original, 'truth', plot_settings, legends, 'o', ...
                                         truth.doriginal, 'truth')
            end
        end
        
        
        function add_individual_plots(obj, plot_t, plot_data, plot_name, plot_settings, legends, marker, dplot_data, dplot_name)
            if nargin <= 7, dplot_data = []; dplot_name = ''; end
            states = sort(plot_settings.States);
            n_cells = min(obj.data.N, plot_settings.MaxCells);
            
            if strcmp(marker, ':'), lw = 1.5; else, lw = 1; end
            
            col_scheme = num2cell(parula(n_cells), 2);
            for k = 1:length(states)
                nexttile(k);
                h = plot(plot_t, reshape(plot_data(:, states(k), 1:n_cells), ...
                                              [], n_cells), marker, 'LineWidth', lw);
                set(h, {'color'}, col_scheme);
                
                if ismember(states(k), obj.data.observed) && states(k) == min(intersect(states, obj.data.observed))
                    add_legendentry(legends(1), h(1), plot_name);
                end
                
                if any(obj.method == "GMGTS") && ~isempty(dplot_data)
                    nexttile(length(states) + k)
                    h = plot(plot_t, reshape(dplot_data(:, states(k), 1:n_cells), ...
                                                  [], n_cells), marker, 'LineWidth', lw);
                    set(h, {'color'}, col_scheme);

                    if ismember(states(k), obj.data.observed) && states(k) == min(intersect(states, obj.data.observed))
                        add_legendentry(legends(2), h(1), dplot_name);
                    end
                end
            end
        end
        
        
        function plot_population_fit(obj, truth, plot_settings)
            states = sort(plot_settings.States);
            n_cells = min(obj.data.N, floor(plot_settings.MaxCells));
            
            screen = get(0, 'ScreenSize');
            figure('position', [55, 25, min(1740, screen(3)-180), min(900, screen(4)-180)])
            tl = tiledlayout(1, length(states));
            title(tl, sprintf('Fitted population trajectories (%d/%d shown)', n_cells, obj.data.N))
            
            for k = 1:length(states)
                nexttile(k)
                title(obj.system.states{states(k)})
                hold on
                if ismember(states(k), obj.data.observed)
                    obs_ind = find(states(k)==obj.data.observed);
                    h = plot(obj.data.t, reshape(obj.data.traces(:, obs_ind, 1:n_cells), ...
                                                  [], n_cells), '-');
                    set(h, {'color'}, num2cell(parula(n_cells), 2));                          
                end
                xlabel('time')
                ylabel('measurement')
                if ismember(states(k), obj.data.observed) && states(k) == min(intersect(states, obj.data.observed))
                    legends = legend(h(1), 'data', 'AutoUpdate', 'off');
                end
            end
            
            mean_label_GMGTS = 'prediction';
            mean_label_GTS = 'prediction';
            CI_label_GMGTS = '90\% CI';
            CI_label_GTS = '90\% CI';
            if numel(obj.method) == 2
                mean_label_GMGTS = [mean_label_GMGTS ' (GMGTS)'];
                mean_label_GTS = [mean_label_GTS ' (GTS)'];
                CI_label_GMGTS = [CI_label_GMGTS ' (GMGTS)'];
                CI_label_GTS = [CI_label_GTS ' (GTS)'];
            end
            
            if ismember("GMGTS", obj.method)
                obj.add_population_plot(obj.results_GMGTS.t_fine, obj.results_GMGTS.population, mean_label_GMGTS, ...
                                        obj.results_GMGTS.fitted_ss, CI_label_GMGTS, ...
                                        plot_settings, legends, ':')
            end
            if ismember("GTS", obj.method)
                obj.add_population_plot(obj.results_GTS.t_fine, obj.results_GTS.population, mean_label_GTS, ...
                                        obj.results_GTS.fitted_ss, CI_label_GTS, ...
                                        plot_settings, legends, '-.')
            end
            if ~isempty(fieldnames(truth))
%                 population = obj.system.integrate(truth.b, obj.data);
                obj.add_population_plot(obj.data.t, truth.moriginal, 'truth', ...
                                        truth.original, '90\% CI (true)', ...
                                        plot_settings, legends, '--')
            end
        end
        
        
        function add_population_plot(obj, plot_t, mean_fitted, mean_label, sample_fitted, ...
                                     sample_label, plot_settings, legends, marker)
            states = sort(plot_settings.States);
            CI = quantile(sample_fitted, [.05 .95], 3);
            switch marker
                case ':', col = [0, 0.4470, 0.7410];
                case '--', col = [0.8500, 0.3250, 0.0980];
                case '-.', col = [0.4940, 0.1840, 0.5560];
            end
            
            for k = 1:length(states)
                nexttile(k)
                h1 = plot(plot_t, mean_fitted(:, states(k)), ...
                          marker, 'LineWidth', 2, 'Color', col);
                plot(plot_t, CI(:, states(k), 1), marker, 'LineWidth', 1.5, 'Color', col);
                plot(plot_t, CI(:, states(k), 2), marker, 'LineWidth', 1.5, 'Color', col);
                h2 = patch([plot_t flip(plot_t)], ...
                           [CI(:, states(k), 1)' flip(CI(:, states(k), 2))'], ...
                           col, 'FaceAlpha', .1, 'LineStyle', 'none');
                if ismember(states(k), obj.data.observed) && states(k) == min(intersect(states, obj.data.observed))
                    add_legendentry(legends, [h1 h2], {mean_label sample_label});
                end
            end
        end
    
    
        function plot_mixed_effects(obj, truth, plot_settings)
            [params, ~, ind] = intersect(sort(plot_settings.Parameters), 1:obj.system.P);
            n_params = length(params);
            
            screen = get(0, 'ScreenSize');
            figure('position', [70, 10, min(1740, screen(3)-180), min(900, screen(4)-180)])
            tl = tiledlayout(length(params), length(params));
            title(tl, 'Random effects marginal distributions')
            
            legends = [];
            parameter_names = obj.system.parameters(~ismember(obj.system.parameters, obj.system.fixed.names));
            for p = 1:n_params
                for q = 1:p
                    nexttile(q + (p-1)*n_params)
                    hold on
                    if q == 1
                        ylabel(parameter_names(ind(p)))
                        if p == 1 || p == n_params
                            legends = [legends legend('AutoUpdate', 'off')]; %#ok<AGROW>
                        end
                    end
                    if p == n_params
                        xlabel(parameter_names(ind(q)))
                    end
                end
            end
            
            labels_ondiag_GMGTS = {'$\hat{\textrm{\boldmath$\beta$}}_i$', 'N$(\hat{\mathbf b}, \hat{\mathbf D})$'};
            labels_ondiag_GTS = {'$\hat{\textrm{\boldmath$\beta$}}_i$', 'N$(\hat{\mathbf b}, \hat{\mathbf D})$'};
            labels_offdiag_GMGTS = {'$\textrm{\boldmath$\beta$}_i$', '$\mathbf b$', '$\mathbf D$'};
            labels_offdiag_GTS = {'$\textrm{\boldmath$\beta$}_i$', '$\mathbf b$', '$\mathbf D$'};
            labels_ondiag_true = strcat(labels_ondiag_GMGTS, ' (true)');
            labels_offdiag_true = strcat(labels_offdiag_GMGTS, ' (true)');
            if numel(obj.method) == 2
                labels_ondiag_GMGTS = strcat(labels_ondiag_GMGTS, ' (GMGTS)');
                labels_ondiag_GTS = strcat(labels_ondiag_GTS, ' (GTS)');
                labels_offdiag_GMGTS = strcat(labels_offdiag_GMGTS, ' (GMGTS)');
                labels_offdiag_GTS = strcat(labels_offdiag_GTS, ' (GTS)');
            end
            
            if ismember("GMGTS", obj.method)
                obj.add_mixed_effects(obj.results_GMGTS.beta_fs(:, ind), ...
                                      obj.results_GMGTS.b_est(ind), ...
                                      obj.results_GMGTS.D_est(ind, ind), ...
                                      obj.results_GMGTS.lognormal, ...
                                      labels_ondiag_GMGTS, labels_offdiag_GMGTS, ...
                                      plot_settings, legends, '^')
            end
            if ismember("GTS", obj.method)
                obj.add_mixed_effects(obj.results_GTS.beta_fs(:, ind), ...
                                      obj.results_GTS.b_est(ind), ...
                                      obj.results_GTS.D_est(ind, ind), ...
                                      obj.results_GTS.lognormal, ...
                                      labels_ondiag_GTS, labels_offdiag_GTS, ...
                                      plot_settings, legends, 'square')
            end
            if ~isempty(fieldnames(truth))
                obj.add_mixed_effects(truth.beta(:, ind), ...
                                      truth.b(ind), ...
                                      truth.D(ind, ind), ...
                                      truth.lognormal, ...
                                      labels_ondiag_true, labels_offdiag_true, ...
                                      plot_settings, legends, 'pentagram')
            end
        end
    

        function add_mixed_effects(obj, beta, b, D, lognormal, labels_ondiag, labels_offdiag, plot_settings, legends, marker)
            params = intersect(sort(plot_settings.Parameters), 1:obj.system.P);
            
            switch marker
                case '^', col = [0, 0.4470, 0.7410]; center = ':v'; lw = 1.2;
                case 'pentagram', col = [0.8500, 0.3250, 0.0980]; center = '--hexagram'; lw = 1;
                case 'square', col = [0.4940, 0.1840, 0.5560]; center = '-.diamond'; lw = 1;
            end
            
            n_params = length(params);
            for p = 1:n_params
                for q = 1:p
                    nexttile(q + (p-1)*n_params)
                    if q == p
                        h1 = histogram(beta(:, p), 20, 'Normalization', 'pdf', 'FaceColor', col, 'FaceAlpha', .3);
                        grid = linspace(min(beta(:, p)), max(beta(:, p)), 30);
                        if lognormal
                            pdf = lognpdf(grid, b(p), sqrt(D(p, p)));
                        else
                            pdf = normpdf(grid, b(p), sqrt(D(p, p)));
                        end
                        h2 = plot(grid, pdf, 'Color', col, 'LineWidth', 1.5);
                        patch([grid flip(grid)], [pdf 0*pdf], col, 'FaceAlpha', .4);
                        if p == 1
                            add_legendentry(legends(1), [h1 h2], labels_ondiag);
                        end
                    else
%                         disp('')
                        h1 = scatter(beta(:, q), beta(:, p), [], col, marker, 'MarkerEdgeAlpha', .3, 'LineWidth', lw);
                        if lognormal
                            h2 = Estimator.mvlncontour(b([q p]), D([q p], [q p]), center, 'Color', col, 'LineWidth', lw);
                        else
                            h2 = Estimator.mvncontour(b([q p]), D([q p], [q p]), center, 'Color', col, 'LineWidth', lw);
                        end
                        if p == n_params && q == 1
                            add_legendentry(legends(2), [h1 h2], labels_offdiag);
                        end
                    end
                end
            end
        end
    end
    
        
    methods (Static)
        function h = mvncontour(m, S, marker, levels, varargin)
            manual = true;
            if ~isnumeric(levels)
                varargin = [{levels} varargin];
                levels = [.68 .95];
                manual = false;
            end
            
            ws = warning('off', 'MATLAB:nearlySingularMatrix');
            
            h1 = plot(m(1), m(2), marker, varargin{:});
%             h1 = plot(m(1), m(2), marker, 'Color', [1 0 1], varargin{:});
            
            S = nearestSPD(S);
            x = linspace(m(1) - 3*sqrt(S(1, 1)), m(1) + 3*sqrt(S(1, 1)), 40);
            y = linspace(m(2) - 3*sqrt(S(2, 2)), m(2) + 3*sqrt(S(2, 2)), 40);

            p = zeros(length(x), length(y));
            for i = 1:length(x)
                for j = 1:length(y)
                    p(i, j) = mvnpdf([x(i) y(j)], m, S);
                end
            end

            S_inv = pinv(S);
            quantiles = sqrt(chi2inv(levels, 2) / S_inv(1, 1));
            levels = mvnpdf(m+[quantiles' zeros(length(levels), 1)], m, S);
%             if max(levels) > 1e8, levels = min(levels, 1e16) / max(levels) * 1e8; end
            if any(~isfinite(levels))
                disp(1)
            end
            if manual
                [~, h2] = contour(kron(x', ones(1, length(y))), kron(y, ones(length(x), 1)), p, marker);
%                 [~, h2] = contour(kron(x', ones(1, length(y))), kron(y, ones(length(x), 1)), p, marker, varargin{:});
            else
                [~, h2] = contour(kron(x', ones(1, length(y))), kron(y, ones(length(x), 1)), p, marker, varargin{:});
            end
            h2.LevelList = levels;
            
            h = [h1 h2];
            
            warning(ws)
        end
        
        
        function h = mvlncontour(m, S, marker, levels, varargin)
            manual = true;
            if ~isnumeric(levels)
                varargin = [{levels} varargin];
                levels = [.68 .95];
%                 levels = (normcdf(.25:.25:2)-.5)*2;
                manual = false;
            end
            
            ws = warning('off', 'MATLAB:nearlySingularMatrix');
            
            h1 = plot(exp(m(1)), exp(m(2)), marker, varargin{:});
%             h1 = plot(m(1), m(2), marker, 'Color', [1 0 1], varargin{:});
            
            S = nearestSPD(S);
            x = linspace(exp(m(1) - 3*sqrt(S(1, 1))), exp(m(1) + 3*sqrt(S(1, 1))), 40);
            y = linspace(exp(m(2) - 3*sqrt(S(2, 2))), exp(m(2) + 3*sqrt(S(2, 2))), 40);

            p = zeros(length(x), length(y));
            for i = 1:length(x)
                for j = 1:length(y)
                    p(i, j) = Estimator.mvlnpdf([x(i) y(j)], m, S);
                end
            end

            S_inv = pinv(S);
            quantiles = sqrt(chi2inv(levels, 2) / S_inv(1, 1));
            levels_transf = Estimator.mvlnpdf(exp(m+[quantiles' zeros(length(levels), 1)]), m, S);
%             if max(levels) > 1e8, levels = min(levels, 1e16) / max(levels) * 1e8; end
            if any(~isreal(p), 'all')
                disp(1)
            end
            if manual
                [~, h2] = contour(kron(x', ones(1, length(y))), kron(y, ones(length(x), 1)), p, marker);
%                 [~, h2] = contour(kron(x', ones(1, length(y))), kron(y, ones(length(x), 1)), p, marker, varargin{:});
            else
                [~, h2] = contour(kron(x', ones(1, length(y))), kron(y, ones(length(x), 1)), p, marker, varargin{:});
            end
            h2.LevelList = levels_transf;
            
            h = [h1 h2];
            
            warning(ws)
        end
        
        
        function p = mvlnpdf(x, m, S)
            [n, d] = size(x);
            p = zeros(size(x, 1), 1);
            for i = 1:n
                p(i) = (2*pi)^(-d/2) * det(S)^-1 * prod(x(i, :)) ...
                                     * exp(-.5 * (log(x(i, :))-m) * tryinv(S) * (log(x(i, :))-m)');
            end
        end
    end
end