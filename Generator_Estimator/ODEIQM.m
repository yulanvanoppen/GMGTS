classdef ODEIQM < handle
    properties (Access = private)
        df_cell                                                             % raw jacobian handle
        g_cell                                                              % rhs: f(x; p) = g(x) p
        dg_cell                                                             % d/dx of rhs
        const_handle                                                               % constant part of rhs
        dconst_cell                                                            % constant part gradient
        dfdx_handle
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
        fixed                                                               % fixed parameter names and values
        fixed_indices                                                       % fixed parameter indices
        x0                                                                  % model initial conditions
        k0                                                                  % model initial parameters
    end
    
    methods
        function obj = ODEIQM(model_file, varargin)                     % Constructor
            parser = inputParser();
            addRequired(parser, 'model_file', @(x) (ischar(x) || isstring(x)) ...
                                                   && endsWith(x, '.txt') );
            addParameter(parser, 'FixedParameters', string([]), @isstring);
            addParameter(parser, 'FixedValues', [], @isnumeric);
            parse(parser, model_file, varargin{:});
            
            assert(isempty(parser.Results.FixedValues) || ...
                   numel(parser.Results.FixedParameters) == numel(parser.Results.FixedValues))
            
            obj.fixed = struct('names', parser.Results.FixedParameters, ...
                               'values', parser.Results.FixedValues);
            obj.process(model_file);
        end
        
        
        function traces = integrate(obj, parameters, settings, t, tol)       % Numerically integrate ODE system
            if nargin < 4, t = settings.t; end
            if nargin < 5, tol = 1e-6; end
            options = struct('reltol', tol);
            
            nseries = size(parameters, 1);
            
            if size(settings.init, 1) == 1
                settings.init = repmat(settings.init, nseries, 1);       % replicate if necessary
            end
            
            full_parameters = zeros(1, length(obj.parameters));
            full_parameters(obj.fixed_indices) = obj.fixed.values;
            traces = zeros(length(t), obj.K, nseries);             % initialize trajectories
            for cell = 1:nseries
                full_parameters(~ismember(1:length(obj.parameters), obj.fixed_indices)) = parameters(cell, :);
                out = obj.MEXf(t, settings.init(cell, :), ...               % integrate using MEXmodel
                               full_parameters, options);
                traces(:, :, cell) = out.statevalues;                       % save integrated states
            end
        end
        
        
        function f = rhs(obj, states, times, parameters)                       % RHS constructed from g()                          
            if nargin < 4, parameters = times; times = zeros(1, size(states, 1)); end
            par_reshaped = reshape(parameters', obj.P, 1, []);              % organize into pages
            f_reshaped = pagemtimes(obj.g(states, times), par_reshaped);           % multiply matrices page-wise
            f = reshape(f_reshaped, size(states)) + obj.const(states, times);      % return right-hand side
        end
        
        
        function S = sensitivity(obj, parameters, settings, t, tol)
            if nargin < 4, t = settings.t; end
            if nargin < 5, tol = 1e-6; end
            options = struct('reltol', tol);
            
            nseries = size(parameters, 1);
            
            if size(settings.init, 1) == 1
                settings.init = repmat(settings.init, nseries, 1);       % replicate if necessary
            end
            
            full_parameters = zeros(1, length(obj.parameters));
            full_parameters(obj.fixed_indices) = obj.fixed.values;
            
            S = zeros(length(t), obj.K, obj.P, nseries);
            for cell = 1:nseries
                full_parameters(~ismember(1:length(obj.parameters), obj.fixed_indices)) = parameters(cell, :);
                out = IQMsensitivity(obj.MEX, t, cellstr(obj.parameters_variable'), [], options, full_parameters, settings.init(cell, :));
                S(:, :, :, cell) = reshape([out.paramtrajectories.states{:}], [], obj.K, obj.P);
            end
        end
        
        
        function result = df(obj, states, times, parameters)
            assert(size(states, 2) == obj.K)
            assert(size(states, 1) == numel(times))
            assert(size(parameters, 2) == obj.P)
            result = cell(obj.K, 1);
            for k = 1:obj.K
                result{k} = obj.df_cell{k}(states, times, parameters);
            end
            result = vertcat(result{:});
        end
        
        
        function result = dfdx(obj, states, dstates, times, parameters)
            assert(size(states, 2) == obj.K)
            assert(size(dstates, 2) == obj.K)
            assert(size(states, 1) == numel(times))
            assert(size(parameters, 2) == obj.P)
            result = obj.dfdx_handle(states, dstates, times, parameters);
%             result = reshape(permute(result, [2 1 3]), obj.P, [], obj.K, size(result, 3));
%             result = permute(result, [2 3 1 4]);
        end
        
        function result = g(obj, states, times)
            assert(size(states, 2) == obj.K)
            assert(size(states, 1) == numel(times))
            result = cell(obj.K, 1);
            for k = 1:obj.K
                result{k} = obj.g_cell{k}(states, times);
            end
            result = vertcat(result{:});
        end
        
        
        function result = dg(obj, states, times)
            assert(size(states, 2) == obj.K)
            assert(size(states, 1) == numel(times))
            result = cell(obj.K, 1);
            for k = 1:obj.K
                result{k} = obj.dg_cell{k}(states, times);
            end
            result = permute(cat(5, result{:}), [1 2 3 5 4]);
        end
        
        
        function result = const(obj, states, times)
            assert(size(states, 2) == obj.K)
            assert(size(states, 1) == numel(times))
            result = obj.const_handle(states, times);
        end
        
        
        function result = dconst(obj, states, times)
            assert(size(states, 2) == obj.K)
            assert(size(states, 1) == numel(times))
            result = cell(obj.K, 1);
            for k = 1:obj.K
                result{k} = obj.dconst_cell{k}(states, times);
            end
            result = horzcat(result{:});
        end
        
        
        
        
        function process(obj, model_file)
            tic
            
            
            %% process model file
            contents = fileread(mode_file);                                 % extract system name manually
            token = regexp(contents, '\*\*\*\*\*\*\*\*\*\* MODEL NAME\s*([^\n]+)\s*\*\*\*\*\*\*\*\*\*\* MODEL NOTES', 'tokens');
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
            
            %% setup symbolic function
            obj.K = length(obj.states);
            t = sym('t', 'real');
            x = sym('x', [obj.K, 1], 'real');
            dx = sym('dx', [obj.K, 1], 'real');
%             f = subs(str2sym(f), cellstr(obj.states'), x');
            f = subs(str2sym(f), [cellstr(obj.states') {'time'}], [x' t]);
            
            intersection = ismember(obj.fixed.names, obj.parameters);
            obj.fixed.names = obj.fixed.names(intersection);
            if ~isempty(obj.fixed.values)
                obj.fixed.values = obj.fixed.values(intersection);
            else
                [~, ind] = ismember(obj.fixed.names, obj.parameters);
                obj.fixed.values = obj.k0(ind)';
            end
            f = subs(f, cellstr(obj.fixed.names), obj.fixed.values);
            
            variable_elements = ~ismember(obj.parameters, obj.fixed.names);
            obj.parameters_variable = obj.parameters(variable_elements);
            obj.k0 = obj.k0(variable_elements);
            obj.P = sum(variable_elements);
            beta = sym('beta', [obj.P, 1], 'real');
            f = subs(f, cellstr(obj.parameters_variable'), beta');
            
            %% linear part
            const_beta = subs(f, beta, zeros(obj.P, 1))
            
            g_symb = sym('g_symb', [obj.K, obj.P], 'real');
            for p = 1:obj.P
                g_symb(:, p) = subs(f, beta, 1*(p == 1:obj.P)') - const_beta;
            end
            for n = 1:numel(g_symb)
                if hasSymType(g_symb(n), 'constant') && g_symb(n) == 0
                    g_symb(n) = 1e-16*x(1);
                elseif ~hasSymType(g_symb(n), 'variable') && g_symb(n) ~= 0
                    g_symb(n) = subs(g_symb(n)) + 1e-16*x(1);
                end
            end
            
            for l = 1:obj.K
                g = matlabFunction(g_symb(l, :), 'Vars', [x; t]);
                g = @(states, times) g(states{:}, times{:});
                obj.g_cell{l} = @(states, times) g(mat2cell(states, size(states, 1), ...
                                                            ones(1, obj.K), size(states, 3)), ...
                                                   mat2cell(repmat(reshape(times, [], 1), 1, 1, size(states, 3)), ...
                                                            numel(times), 1, size(states, 3)));
                
                dg = jacobian(reshape(g_symb, [], 1), x(l));
                dg_vectorized = dg;
                for n = 1:numel(dg_vectorized)
                    if dg_vectorized(n) == 0
                        dg_vectorized(n) = 1e-16*x(1);
                    elseif ~hasSymType(dg_vectorized(n), 'variable') && dg_vectorized(n) ~= 0
                        dg_vectorized(n) = subs(dg_vectorized(n)) + 1e-16*x(1);
                    end
                end
                dg = matlabFunction(dg_vectorized, 'Vars', [x; t]);
                dg = @(states, times) dg(states{:}, times{:});
                dg = @(states, times) dg(mat2cell(permute(reshape(permute(states, [2 1 3]), ...
                                                           obj.K, 1, size(states, 1), ...
                                                           size(states, 3)), [2 1 3 4]), ...
                                           1, ones(1, obj.K), size(states, 1), size(states, 3)), ...
                                         mat2cell(repmat(reshape(times, 1, 1, []), 1, 1, 1, size(states, 3)), ...
                                                            1, 1, numel(times), size(states, 3)));
                obj.dg_cell{l} = @(states, times) reshape(dg(states, times), obj.K, obj.P, ...
                                                          size(states, 1), size(states, 3));
                                                      
                dconst_beta = jacobian(const_beta, x(l));
                dconst_vectorized = dconst_beta;
                for n = 1:numel(dconst_vectorized)
                    if dconst_vectorized(n) == 0
                        dconst_vectorized(n) = 1e-16*x(1);
                    elseif ~hasSymType(dconst_vectorized(n), 'variable') && dconst_vectorized(n) ~= 0
                        dconst_vectorized(n) = subs(dconst_vectorized(n)) + 1e-16*x(1);
                    end
                end
                dconst_beta = matlabFunction(dconst_vectorized, 'Vars', [x; t]);
                dconst_beta = @(states, times) dconst_beta(states{:}, times{:});
                obj.dconst_cell{l} = @(states, times) dconst_beta(mat2cell(states, size(states, 1), ...
                                                                           ones(1, obj.K), size(states, 3)), ...
                                                                  mat2cell(repmat(reshape(times, [], 1), 1, 1, size(states, 3)), ...
                                                                           numel(times), 1, size(states, 3)));
            end
            
            
            
            for n = 1:numel(const_beta)
                if hasSymType(const_beta(n), 'constant') && const_beta(n) == 0
                    const_beta(n) = 1e-16*x(1);
                elseif ~hasSymType(const_beta(n), 'variable') && const_beta(n) ~= 0
                    const_beta(n) = subs(const_beta(n)) + 1e-16*x(1);
                end
            end
            
            const_beta = matlabFunction(permute(const_beta, [2 1]), 'Vars', [x; t]);
            const_beta = @(states, times) const_beta(states{:}, times{:});
            obj.const_handle = @(states, times) const_beta(mat2cell(states, size(states, 1), ...
                                                                    ones(1, obj.K), size(states, 3)), ...
                                                           mat2cell(repmat(reshape(times, [], 1), 1, 1, size(states, 3)), ...
                                                                    numel(times), 1, size(states, 3)));
            
            %% RHS derivatives
            df_symb = jacobian(f, x)
            for n = 1:numel(df_symb)
                if hasSymType(df_symb(n), 'constant') && df_symb(n) == 0
                    df_symb(n) = 1e-16*x(1);
                elseif ~hasSymType(df_symb(n), 'variable') && df_symb(n) ~= 0
                    df_symb(n) = subs(df_symb(n)) + 1e-16*x(1);
                end
            end
            for l = 1:obj.K
                df = matlabFunction(df_symb(l, :), 'Vars', [x; t; beta]);
                df = @(states, times, parameters) df(states{:}, times{:}, parameters{:});
                obj.df_cell{l} = @(states, times, parameters) ...
                    df(mat2cell(states, size(states, 1), ...
                                ones(1, obj.K), size(states, 3)), ...
                       mat2cell(repmat(reshape(times, [], 1), 1, 1, size(states, 3)), ...
                                numel(times), 1, size(states, 3)), ...
                       mat2cell(repmat(permute(parameters, [3 2 1]), size(states, 1), 1, 1), ...
                                size(states, 1), ones(1, obj.P), size(parameters, 1)));
            end
            
            dfxdot = (jacobian(f, x) * dx)';
            dfxdot_vectorized = dfxdot;
            for n = 1:numel(dfxdot_vectorized)
                if dfxdot_vectorized(n) == 0
                    dfxdot_vectorized(n) = 1e-16*dx(1);
                elseif ~hasSymType(dfxdot_vectorized(n), 'variable') && dfxdot_vectorized(n) ~= 0
                    dfxdot_vectorized(n) = subs(dfxdot_vectorized(n)) + 1e-16*dx(1);
                end
            end
            dfxdot = matlabFunction(dfxdot_vectorized, 'Vars', [x; dx; t; beta]);
            dfxdot = @(states, dstates, times, parameters) dfxdot(states{:}, dstates{:}, times{:}, parameters{:});
            obj.dfdx_handle = @(states, dstates, times, parameters) ...
                permute(dfxdot(mat2cell(repmat(reshape(states, size(states, 1), size(states, 2), 1, size(states, 3)), ...
                                               1, 1, size(dstates, 3), 1), ...
                                        size(states, 1), ones(1, obj.K), size(dstates, 3), size(states, 3)), ...
                               mat2cell(dstates, size(dstates, 1), ones(1, obj.K), size(dstates, 3), size(dstates, 4)), ...
                               mat2cell(repmat(reshape(times, [], 1), 1, 1, size(dstates, 3), size(states, 3)), ...
                                        numel(times), 1, size(dstates, 3), size(states, 3)), ...
                               mat2cell(repmat(permute(parameters, [3 2 4 1]), size(states, 1), 1, size(dstates, 3), 1), ...
                                        size(states, 1), ones(1, obj.P), size(dstates, 3), size(parameters, 1))), [1 2 3 4]);
            
            toc_converting = toc
        end
    end
end