classdef Estimator < handle
    properties (Access = public)
        data                                                                % observed data struct
        system                                                              % ODE system object
        stages                                                              % stages to execute
        method                                                              % method(s) to conduct inference
        knots
        linear
        lambda
        lb
        ub
        prior
        state_weights
        
        GMGTS_settings                                                      % GMGTS settings struct
        GMGTS_smoother                                                      % GMGTS smoothing object
        GMGTS_first_stage                                                   % GMGTS first stage optimization object
        GMGTS_second_stage                                                  % GMGTS second stage inference object
        results_GMGTS                                                       % collect GMGTS results
        
        GTS_settings                                                        % GTS settings struct
        GTS_first_stage                                                     % GTS first stage optimization object
        GTS_second_stage                                                    % GTS second stage inference object
        results_GTS                                                         % collect GTS results
    end
    
    methods
        %% Constructor -----------------------------------------------------
        function obj = Estimator(data, system, stages, varargin)
            default_Stages = 2;
            default_Methods = "GMGTS";
            
            initial = system.k0';
            default_Knots = linspace(data.t(1), data.t(end), round((data.T-1)*.75)+1);
            default_Knots = default_Knots(2:end-1);
            default_Linear = [0 10000];
            default_LB = .25 * initial;
            default_UB = 4 .* initial + .0001 * mean(initial);
            default_Prior = struct('mean', 0, 'prec', 0);
            default_StateWeights = ones(1, system.K);
            
            parser = inputParser;
            addRequired(parser, 'data', @isstruct);
            addRequired(parser, 'system', @(x) isa(x, 'ODEIQM'));
            addParameter(parser, 'Stages', default_Stages, @(x) ismember(x, 0:2));
            addParameter(parser, 'Methods', default_Methods, ...
                         @(x) all(ismember(string(x), ["GMGTS" "GTS"])));
            addParameter(parser, 'Knots', default_Knots, @(x) all(data.t(1) <= x & x <= data.t(end)));
            addParameter(parser, 'Linear', default_Linear, @(x) numel(x) == 2 && x(1) < x(2));
            addParameter(parser, 'Lambda', [], @(x) all(x > 0))
            addParameter(parser, 'LB', default_LB, @(x) all(x < initial));
            addParameter(parser, 'UB', default_UB, @(x) all(x > initial));
            addParameter(parser, 'LB_States', default_LB, @(x) all(x < initial));
            addParameter(parser, 'UB_States', default_UB, @(x) all(x > initial));
            addParameter(parser, 'Prior', default_Prior, @(x) isfield(x, 'mean') && numel(x.mean) == system.P ...
                                                           && (isfield(x, 'prec') && all(size(x.prec) == system.P) ...
                                                               && issymmetric(x.prec) && all(eig(x.prec) >= 0) ...
                                                            || isfield(x, 'sd') && numel(x.sd) == system.P ...
                                                               && all(x.sd > 0) ...
                                                            || isfield(x, 'cv') && numel(x.cv) == system.P ...
                                                               && all(x.cv > 0)));
            addParameter(parser, 'StateWeights', default_StateWeights, @(x) numel(x) == system.K && all(x > 0));
            parse(parser, data, system, stages, varargin{:});
            
            obj.data = parser.Results.data;
            obj.system = parser.Results.system;
            obj.stages = parser.Results.Stages;
            obj.method = string(parser.Results.Methods);
            obj.knots = ( sort(unique([data.t(1) parser.Results.Knots data.t(end)])) - data.t(1) ) / range(data.t);
            obj.linear = parser.Results.Linear;
            obj.lambda = parser.Results.Lambda;
            obj.lb = parser.Results.LB;
            obj.ub = parser.Results.UB;
            obj.state_weights = parser.Results.StateWeights;

            obj.prior = parser.Results.Prior;
            obj.prior.mean = reshape(obj.prior.mean, [], 1);
            if isfield(obj.prior, 'cv')
                obj.prior.cv = reshape(obj.prior.cv, [], 1);
                obj.prior.mean(obj.prior.mean == 0 & isinf(obj.prior.cv)) = 1;
                obj.prior.sd = obj.prior.cv .* obj.prior.mean;
            end
            if isfield(obj.prior, 'sd')
                obj.prior.sd = reshape(obj.prior.sd, [], 1);
                warning('off','MATLAB:singularMatrix')
                obj.prior.prec = inv(diag(obj.prior.sd.^2));
                warning('on','MATLAB:singularMatrix')
            end
            if ~isfield(obj.prior, 'mult'), obj.prior.mult = 1; end
            
            obj.state_weights = reshape(obj.state_weights, 1, []) / sum(obj.state_weights, 'all');
            
            if ismember("GMGTS", obj.method), obj.constructor_GMGTS; end
            if ismember("GTS", obj.method), obj.constructor_GTS; end
        end
        
        
        function constructor_GMGTS(obj)
%             grid = linspace(obj.data.t(1), 0.75*obj.data.t(end), 11);
            grid = obj.data.t(1) + (0:.1:1) * range(obj.data.t);
%             grid = obj.data.t(1) + (0:.1:1) * range(obj.data.t)/2;
            weights = ones(length(grid), obj.system.K);
            weights([1 end], :) = 0;
%             optindices = [1:11 13:2:21 25:4:40];
            optindices = 1:length(obj.data.t);
            
%             disp(obj.linear)

%             interactive = [true true];
%             interactive = [false true];
%             interactive = [true false];
            interactive = [false false];
            obj.GMGTS_settings.sm = struct('order', 4, 'knots', obj.knots, 'initial', obj.system.k0', ...
                                           'linear', obj.linear, 'lambda', obj.lambda, 'grid', grid, 'interactive', interactive(1));
            obj.GMGTS_settings.fs = struct('grid', grid, 'weights', weights, 'initial', obj.system.k0', ...
                                           'nrep', 5, 'tol', 2e-3, 'nstart', 10, 'perturbation', 0, ...
                                           'lb', obj.lb, 'ub', obj.ub, 'optindices', optindices, 'prior', obj.prior, ...
                                           'state_weights', obj.state_weights, 'interactive', interactive(2));
            obj.GMGTS_settings.ss = struct('weights', weights, 'nrep', 10, ...
                                           'tol', 1e-2, 'prior', obj.prior);
        end
        
        
        function constructor_GTS(obj)
%             optindices = 1:16;
            optindices = 1:length(obj.data.t);
            obj.GTS_settings.fs = struct('initial', obj.system.k0', 'lb', obj.lb, 'ub', obj.ub, 'optindices', optindices, ...
                                         'nrep', 5, 'tol', 2e-3, 'nstart', 10, 'perturbation', 0, 'prior', obj.prior);
            obj.GTS_settings.ss = struct('nrep', 10, 'tol', 1e-2, 'prior', obj.prior);
        end
        
        
        %% Estimation ------------------------------------------------------
        function estimate(obj, silent)
            if nargin < 2, silent = false; else, silent = true; end
            if ismember("GMGTS", obj.method), obj.estimate_GMGTS(silent); end
            if ismember("GTS", obj.method), obj.estimate_GTS(silent); end
        end
        
        
        function estimate_GMGTS(obj, silent)
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
        end
        
        
        function estimate_GTS(obj, silent)
            tic
            obj.GTS_first_stage = FirstStageGTS(obj.data, obj.system, obj.GTS_settings.fs);
            obj.results_GTS = obj.GTS_first_stage.optimize();
            if ~silent, toc_fs_GTS = toc, else, toc_fs_GTS = toc; end
            toc_ss_GTS = 0;

            if obj.stages == 2
                tic
                obj.GTS_second_stage = SecondStageGTS(obj.results_GTS, obj.system, ...
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
                D = obj.results_GMGTS.D_GTS;
                var = @(trace) obj.GMGTS_smoother.theta_fs(1) + obj.GMGTS_smoother.theta_fs(2) * trace.^2;
            else
                b = obj.results_GTS.b_est;
                D = obj.results_GTS.D_GTS;
                var = @(trace) obj.GTS_first_stage.theta_fs(1) + obj.GTS_first_stage.theta_fs(2) * trace.^2;
            end
            
            logfs = @(trace, i) reshape(sum(-1/2 * (obj.data.traces(:, :, i) - trace).^2 ./ var(trace) - 1/2 * log(var(trace)), [1 2]), [], 1);
            
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
            default_States = 1:obj.system.K;
            default_Parameters = 1:obj.system.P;
            default_MaxCells = 5;
            
            overlaps = @(x, S) ~isempty(intersect(string(x), string(S)));
            
            addParameter(parser, 'True', struct, @(x) isstruct(x));
            addParameter(parser, 'States', default_States, ...
                         @(x) overlaps(x, obj.data.observed) || ...
                              overlaps(x, obj.system.states));
            addParameter(parser, 'Parameters', default_Parameters, ...
                         @(x) overlaps(x, obj.data.varying) || ...
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
            
            figure('position', [25, 55, 1740, 900])
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
                h2 = plot(linspace(obj.data.t(1), obj.data.t(end), 201), ...
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
                h4 = plot(linspace(obj.data.t(1), obj.data.t(end), 201), ...
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
            
            figure('position', [40, 40, 1740, 900])
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
                        h = plot(linspace(obj.data.t(1), obj.data.t(end), 201), ...
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
            
            figure('position', [55, 25, 1740, 900])
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
                                        obj.results_GMGTS.fitted2, CI_label_GMGTS, ...
                                        plot_settings, legends, ':')
            end
            if ismember("GTS", obj.method)
                obj.add_population_plot(obj.results_GTS.t, obj.results_GTS.population, mean_label_GTS, ...
                                        obj.results_GTS.fitted2, CI_label_GTS, ...
                                        plot_settings, legends, '-.')
            end
            if ~isempty(fieldnames(truth))
                population = obj.system.integrate(truth.b, obj.data);
                obj.add_population_plot(obj.data.t, population, 'truth', ...
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
            [params, ~, ind] = intersect(sort(plot_settings.Parameters), obj.data.varying);
            n_params = length(params);
            
            figure('position', [70, 10, 1740, 900])
            tl = tiledlayout(length(params), length(params));
            title(tl, 'Random effects marginal distributions')
            
            legends = [];
            parameter_names = obj.system.parameters(~ismember(obj.system.parameters, obj.system.fixed.names));
            for p = 1:n_params
                for q = 1:p
                    nexttile(q + (p-1)*n_params)
                    hold on
                    if q == 1
                        ylabel(parameter_names(obj.data.varying(ind(p))))
                        if p == 1 || p == n_params
                            legends = [legends legend('AutoUpdate', 'off')]; %#ok<AGROW>
                        end
                    end
                    if p == n_params
                        xlabel(parameter_names(obj.data.varying(ind(q))))
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
                obj.add_mixed_effects(obj.results_GMGTS.beta_fs(:, obj.data.varying(ind)), ...
                                      obj.results_GMGTS.b_est(obj.data.varying(ind)), ...
                                      obj.results_GMGTS.D_GTS(ind, ind), ...
                                      labels_ondiag_GMGTS, labels_offdiag_GMGTS, ...
                                      plot_settings, legends, '^')
            end
            if ismember("GTS", obj.method)
                obj.add_mixed_effects(obj.results_GTS.beta_fs(:, obj.data.varying(ind)), ...
                                      obj.results_GTS.b_est(obj.data.varying(ind)), ...
                                      obj.results_GTS.D_GTS(ind, ind), ...
                                      labels_ondiag_GTS, labels_offdiag_GTS, ...
                                      plot_settings, legends, 'square')
            end
            if ~isempty(fieldnames(truth))
                obj.add_mixed_effects(truth.beta(:, obj.data.varying(ind)), ...
                                      truth.b(obj.data.varying(ind)), ...
                                      truth.D(obj.data.varying(ind), obj.data.varying(ind)), ...
                                      labels_ondiag_true, labels_offdiag_true, ...
                                      plot_settings, legends, 'pentagram')
            end
        end
    

        function add_mixed_effects(obj, beta, b, D, labels_ondiag, labels_offdiag, plot_settings, legends, marker)
            params = intersect(sort(plot_settings.Parameters), obj.data.varying);
            
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
                        pdf = normpdf(grid, b(p), sqrt(D(p, p)));
                        h2 = plot(grid, pdf, 'Color', col, 'LineWidth', 1.5);
                        patch([grid flip(grid)], [pdf 0*pdf], col, 'FaceAlpha', .4);
                        if p == 1
                            add_legendentry(legends(1), [h1 h2], labels_ondiag);
                        end
                    else
%                         disp('')
                        h1 = scatter(beta(:, q), beta(:, p), [], col, marker, 'MarkerEdgeAlpha', .3, 'LineWidth', lw);
                        h2 = Estimator.mvncontour(b([q p]), D([q p], [q p]), center, 'Color', col, 'LineWidth', lw);
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
    end
end