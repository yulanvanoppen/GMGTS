classdef System < handle
    properties (Access = private)
        df_cell                                                             % raw jacobian handle
        g_cell                                                              % RHS: f(x; p) = g(x)⋅k + h(x)
        dg_cell                                                             % d/dx of RHS
        h_handle                                                            % constant part of RHS
        dh_cell                                                             % constant part gradient
        dfdx_handle                                                         % RHS Jacobian
        model                                                               % IQMtools model
        MEX                                                                 % MEX model name
        MEXf                                                                % MEX function
    end
    
    properties (SetAccess = private)
        name                                                                % system name
        K                                                                   % system dimension
        P                                                                   % number of system parameters
        states                                                              % state names
        parameters                                                          % parameters
        parameters_variable                                                 % parameters excluding fixed
        fixed_indices                                                       % fixed parameter indices
        variable_indices
        x0                                                                  % model initial conditions
        k0                                                                  % model initial parameters
    end
    
    properties (Access = public)
        fixed                                                               % fixed parameter names and values
    end
    
    
    methods
        function obj = System(model_file, varargin)                     % Constructor
            parser = inputParser;
            parser.KeepUnmatched = true;                                         
            addRequired(parser, 'model_file', @(x) (ischar(x) || isstring(x)) ...
                                                   && endsWith(x, '.txt') );
            addParameter(parser, 'FixedParameters', string([]), @isstring); 
            addParameter(parser, 'FixedValues', [], @isnumeric);
            parse(parser, model_file, varargin{:});
            
            assert(isempty(parser.Results.FixedValues) || ...               % equal number of labels and constraints
                   numel(parser.Results.FixedParameters) == numel(parser.Results.FixedValues))
            
            obj.fixed = struct('names', parser.Results.FixedParameters, ... % store IDs and fixed values
                               'values', parser.Results.FixedValues);
            obj.process(model_file);                                        % extract relevant functions using IQMmodel
        end
        
        
        function traces = integrate(obj, parameters, settings, t, tol)  % Numerically integrate ODE system
            if nargin < 4, t = settings.t; end
            if nargin < 5, tol = 1e-6; end
            options = struct('reltol', tol);
            
            nseries = size(parameters, 1);                                  % number of replicates to integrate
            
            if size(settings.init, 1) == 1
                settings.init = repmat(settings.init, nseries, 1);          % replicate if necessary
            end
            
            full_parameters = zeros(1, length(obj.parameters));             % substitute fixed parameter values
            full_parameters(obj.fixed_indices) = obj.fixed.values;
            
            traces = zeros(length(t), obj.K, nseries);
            for cell = 1:nseries                                            % substitute cell-specific parameters
                full_parameters(obj.variable_indices) = parameters(cell, :);
                out = obj.MEXf(t, settings.init(cell, :), ...               % fast integration using MEXmodel
                               full_parameters, options);
                traces(:, :, cell) = out.statevalues;                       % save integrated states
            end
        end
        
        
        function f = rhs(obj, states, times, parameters)                % RHS constructed from g()                          
            if nargin < 4, parameters = times; times = zeros(1, size(states, 1)); end
            par_reshaped = reshape(parameters', obj.P, 1, []);              % organize into pages
            f_reshaped = pagemtimes(obj.g(states, times), par_reshaped);    % multiply matrices page-wise
            f = reshape(f_reshaped, size(states)) + obj.h(states, times);   % return right-hand side
        end
        
        
        function S = sensitivity(obj, parameters, settings, t, tol)     % Solve sensitivity equations
            if nargin < 4, t = settings.t; end
            if nargin < 5, tol = 1e-6; end
            options = struct('reltol', tol);
            
            nseries = size(parameters, 1);                                  % number of replicates to integrate
            
            if size(settings.init, 1) == 1
                settings.init = repmat(settings.init, nseries, 1);          % replicate if necessary
            end
            
            full_parameters = zeros(1, length(obj.parameters));             % substitute fixed values
            full_parameters(obj.fixed_indices) = obj.fixed.values;
            
            S = zeros(length(t), obj.K, obj.P, nseries);
            for cell = 1:nseries                                            % substitute cell-specific parameters
                full_parameters(~ismember(1:length(obj.parameters), obj.fixed_indices)) = parameters(cell, :);
                out = IQMsensitivity(obj.MEX, t, cellstr(obj.parameters_variable'), [], options, full_parameters, settings.init(cell, :));
                S(:, :, :, cell) = reshape([out.paramtrajectories.states{:}], [], obj.K, obj.P);
            end
        end
        
        
        function result = df(obj, states, times, parameters)            % Evaluate RHS Jacobian from saved function handles
            assert(size(states, 2) == obj.K)
            assert(size(states, 1) == numel(times))
            assert(size(parameters, 2) == obj.P)
            result = cell(obj.K, 1);
            for k = 1:obj.K                                                 % collect evaluations per partial derivative
                result{k} = obj.df_cell{k}(states, times, parameters);
            end
            result = vertcat(result{:});                                    % stack vertically in K (TxKxN-dimensional) blocks 
        end
        
        
        function result = g(obj, states, times)                         % Evaluate RHS matrix factor from saved function handles
            assert(size(states, 2) == obj.K)
            assert(size(states, 1) == numel(times))
            result = cell(obj.K, 1);
            for k = 1:obj.K                                                 % collect evaluations per state component
                result{k} = obj.g_cell{k}(states, times);
            end
            result = vertcat(result{:});                                    % stack vertically in K (TxPxN-dimensional) blocks 
        end
        
        
        function result = dg(obj, states, times)                        % Evaluate matrix factor Jacobian from saved function handles
            assert(size(states, 2) == obj.K)
            assert(size(states, 1) == numel(times))
            result = cell(obj.K, 1);
            for k = 1:obj.K                                                 % collect evaluations per partial derivative
                result{k} = obj.dg_cell{k}(states, times);
            end                                                             % concatenate partial derivatives along 4th dimension
            result = permute(cat(5, result{:}), [1 2 3 5 4]);               % in K (KxPxTx1xN-dimensional) blocks
        end
        
        
        function result = h(obj, states, times)                         % Evaluate RHS constant factor from saved function handle
            assert(size(states, 2) == obj.K)
            assert(size(states, 1) == numel(times))
            result = obj.h_handle(states, times);                           % same shape as 'states' input
        end
        
        
        function result = dh(obj, states, times)                        % Evaluate constant factor Jacobian from saved function handles
            assert(size(states, 2) == obj.K)
            assert(size(states, 1) == numel(times))
            result = cell(obj.K, 1);
            for k = 1:obj.K                                                 % collect evaluations per partial derivative
                result{k} = obj.dh_cell{k}(states, times);
            end
            result = horzcat(result{:});                                    % stack vertically in K (TxKxN-dimensional) blocks 
        end
        
        
        function process(obj, model_file)                               % Process model file
            contents = fileread(model_file);                                % regex to extract system name
            token = regexp(contents, '\* MODEL NAME\s*([^\n]+)\s*\*', 'tokens');
            obj.name = strtrim(token{1}{1});
            
            installIQMtools;
            obj.model = IQMmodel(model_file);                               % read model file
            obj.MEX = ['MEX' model_file(6:end-4)];
            IQMmakeMEXmodel(obj.model, obj.MEX);                            % convert to MEXmodel
            obj.MEXf = str2func(obj.MEX);
            
            [obj.states, f, obj.x0] = IQMstates(obj.model);
            [obj.parameters, obj.k0] = IQMparameters(obj.model);
            obj.states = string(obj.states);
            obj.parameters = string(obj.parameters);
            [~, obj.fixed_indices] = ismember(obj.fixed.names, obj.parameters);
            obj.fixed_indices = obj.fixed_indices(obj.fixed_indices ~= 0);
            obj.variable_indices = ~ismember(1:length(obj.parameters), obj.fixed_indices);
            
            [f, t, x, beta] = obj.construct_symbolic(f);
            obj.construct_linear_decomposition(f, t, x, beta);
            obj.construct_rhs_jacobian(f, t, x, beta);
        end
        
        
        function [f, t, x, beta] = construct_symbolic(obj, f)           % Setup symbolic function
            obj.K = length(obj.states);
            t = sym('t', 'real');
            x = sym('x', [obj.K, 1], 'real');
            f = subs(str2sym(f), [cellstr(obj.states') {'time'}], [x' t]);
            
            intersection = ismember(obj.fixed.names, obj.parameters);       % substitute fixed parameters
            obj.fixed.names = obj.fixed.names(intersection);
            if ~isempty(obj.fixed.values)
                obj.fixed.values = obj.fixed.values(intersection);
            else
                [~, ind] = ismember(obj.fixed.names, obj.parameters);
                obj.fixed.values = obj.k0(ind)';
            end
            f = subs(f, cellstr(obj.fixed.names), obj.fixed.values);
            
            variable_elements = ~ismember(obj.parameters, obj.fixed.names); % separate remaining parameters
            obj.parameters_variable = obj.parameters(variable_elements);
            obj.k0 = obj.k0(variable_elements);
            obj.P = sum(variable_elements);
            beta = sym('beta', [obj.P, 1], 'real');
            f = subs(f, cellstr(obj.parameters_variable'), beta');
        end
        
        
        function construct_linear_decomposition(obj, f, t, x, beta)     % Extract RHS decomposition g(x)⋅k+h(x)
            h_symb = subs(f, beta, zeros(obj.P, 1));
            
            g_symb = sym('g_symb', [obj.K, obj.P], 'real');                 % separate h(x) from f(x)
            for p = 1:obj.P
                g_symb(:, p) = subs(f, beta, 1*(p == 1:obj.P)') - h_symb;
            end
            for n = 1:numel(g_symb)
                if hasSymType(g_symb(n), 'constant') && g_symb(n) == 0      % fix for zero entries
                    g_symb(n) = 1e-16*x(1);
                elseif ~hasSymType(g_symb(n), 'variable') && g_symb(n) ~= 0 % fix for constant entries
                    g_symb(n) = subs(g_symb(n)) + 1e-16*x(1);
                end
            end
            
            for l = 1:obj.K                                                 
                g = matlabFunction(g_symb(l, :), 'Vars', [x; t]);           % create vectorized matlabFunction
                g = @(states, times) g(states{:}, times);
                obj.g_cell{l} = @(states, times) g(obj.deal_states(states), ...
                                                   obj.deal_times(times, size(states, 3)));
                
                dg = jacobian(flatten(g_symb), x(l));                       % create vectorized Jacobian as with g
                dg_vectorized = dg;
                for n = 1:numel(dg_vectorized)
                    if dg_vectorized(n) == 0
                        dg_vectorized(n) = 1e-16*x(1);
                    elseif ~hasSymType(dg_vectorized(n), 'variable') && dg_vectorized(n) ~= 0
                        dg_vectorized(n) = subs(dg_vectorized(n)) + 1e-16*x(1);
                    end
                end
                dg = matlabFunction(dg_vectorized, 'Vars', [x; t]);
                dg = @(states, times) dg(states{:}, times);
                dg = @(states, times) dg(obj.deal_states_dg(states), ...
                                         obj.deal_times_dg(times, size(states, 3)));
                obj.dg_cell{l} = @(states, times) reshape(dg(states, times), obj.K, obj.P, ...
                                                          size(states, 1), size(states, 3));
                                                      
                dh_symb = jacobian(h_symb, x(l));                           % create vectorized Jacobian of h similarly
                dh_vectorized = dh_symb;
                for n = 1:numel(dh_vectorized)
                    if dh_vectorized(n) == 0
                        dh_vectorized(n) = 1e-16*x(1);
                    elseif ~hasSymType(dh_vectorized(n), 'variable') && dh_vectorized(n) ~= 0
                        dh_vectorized(n) = subs(dh_vectorized(n)) + 1e-16*x(1);
                    end
                end
                dh_symb = matlabFunction(dh_vectorized, 'Vars', [x; t]);
                dh_symb = @(states, times) dh_symb(states{:}, times);
                obj.dh_cell{l} = @(states, times) dh_symb(obj.deal_states(states), ...
                                                          obj.deal_times(times, size(states, 3)));
            end
            
            for n = 1:numel(h_symb)                                         % fix zero and constant entries
                if hasSymType(h_symb(n), 'constant') && h_symb(n) == 0
                    h_symb(n) = 1e-16*x(1);
                elseif ~hasSymType(h_symb(n), 'variable') && h_symb(n) ~= 0
                    h_symb(n) = subs(h_symb(n)) + 1e-16*x(1);
                end
            end
            
            h_symb = matlabFunction(permute(h_symb, [2 1]), 'Vars', [x; t]);% vectorize analogously
            h_symb = @(states, times) h_symb(states{:}, times);
            obj.h_handle = @(states, times) h_symb(obj.deal_states(states), ...
                                                   obj.deal_times(times, size(states, 3)));
        end
        
        
        function construct_rhs_jacobian(obj, f, t, x, beta)             % RHS state derivatives
            df_symb = jacobian(f, x);                                       % RHS Jacobian
            for n = 1:numel(df_symb)                                        % fix zero and constant entries
                if hasSymType(df_symb(n), 'constant') && df_symb(n) == 0
                    df_symb(n) = 1e-16*x(1);
                elseif ~hasSymType(df_symb(n), 'variable') && df_symb(n) ~= 0
                    df_symb(n) = subs(df_symb(n)) + 1e-16*x(1);
                end
            end
            for l = 1:obj.K
                df = matlabFunction(df_symb(l, :), 'Vars', [x; t; beta]);   % vectorized matlabFunction
                df = @(states, times, parameters) df(states{:}, times, parameters{:});
                obj.df_cell{l} = @(states, times, parameters) df(obj.deal_states(states), ...
                                                                 obj.deal_times(times, size(states, 3)), ...
                                                                 obj.deal_parameters(parameters, size(states, 1)));
            end
        end

        
        function cell_states = deal_states(obj, states)                 % Create cells along state dimension
            cell_states = cell(1, obj.K);
            for k = 1:obj.K, cell_states{k} = states(:, k, :); end
        end

        
        function rep_times = deal_times(~, times, N)                    % Repeat time points for each cell
            rep_times = repmat(flatten(times), 1, 1, N);
        end
        
        
        function cell_parameters = deal_parameters(obj, parameters, T)  % Create cells along parameter component dimension
            rep_parameters = repmat(permute(parameters, [3 2 1]), T, 1, 1); % repeat for each time point
            
            cell_parameters = cell(1, obj.P);
            for p = 1:obj.P, cell_parameters{p} = rep_parameters(:, p, :); end
        end
        
        
        function cell_states = deal_states_dg(obj, states)              % Additionally remove time from first dimension
            [T, ~, N] = size(states);
            transposed = permute(states, [2 1 3]);
            reshaped = reshape(transposed, 1, obj.K, T, N);
            
            cell_states = cell(1, obj.K);
            for k = 1:obj.K, cell_states{k} = reshaped(:, k, :, :); end
        end

        
        function rep_times = deal_times_dg(~, times, N)                 % Repetition along 3rd dimension
            rep_times = repmat(reshape(times, 1, 1, []), 1, 1, 1, N);
        end
    end
end