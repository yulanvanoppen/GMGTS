classdef Smoother < handle
    properties (Access = private)
        T                                                                   % number of time points
        T_fs                                                                % number of optimized time points
        L                                                                   % system dimension
        N                                                                   % number of cells
        penalty_ind
    end
    
    properties (SetAccess = private)
        data                                                                % general data and output container
        system
        settings
        bsplines                                                             % nested object controlling splines from B-splines
        
        delta
        lambda                                                              % penalty multiplier
        B                                                                   % B-spline basis evaluated at t
        dB                                                                  % corresponding first derivative
        ddB                                                                 % corresponding second derivative
        B_fine
        dB_fine
        smoothed                                                            % smoothed data
        dsmoothed                                                           % smoothed derivative
        variances_sm
        variances_fs
        sigma
        tau
    end
    
    methods
        function obj = Smoother(data, system, settings)                 % Constructor
            obj.data = data;                                                % store experiment data
            obj.system = system;
            obj.settings = settings;
            [obj.T, obj.L, obj.N] = size(obj.data.traces);                  % extract trajectory dimensions
            
            obj.bsplines = cell(1, obj.L);
            obj.delta = cell(1, obj.L);
            for state = 1:obj.L                                             % instantiate spline
                obj.bsplines{state} = BSpline(settings.order, settings.knots{state});
                obj.delta{state} = zeros(obj.bsplines{state}.card, obj.N);
            end
            
            if ~isempty(settings.lambda)                                    % initial penalty multiplier
                obj.lambda = reshape(settings.lambda, 1, []) .* ones(1, obj.L);
            else
                obj.lambda = 1e-3 * ones(1, obj.L);
            end
            
            obj.variances_sm = zeros(obj.T, obj.L, obj.N);
            obj.variances_fs = zeros(length(settings.grid), obj.L, obj.N);  % measurement error variances
            obj.sigma = zeros(1, obj.L);
            obj.tau = zeros(1, obj.L);
        end
        
        
        function output = smooth(obj)                                   % Smooth trajectory data
            if obj.settings.interactive
                app = smoothing_app(obj);
                waitfor(app.FinishButton, 'UserData')
                if isvalid(app)
                    obj.data = app.smoother.data;
                    delete(app)
                else
                    obj.optimize()
                end
            else
                obj.optimize()                                              % return extended data struct
            end
            
            obj.data.t_data = obj.data.t;
            obj.data.T_data = obj.data.T;
            obj.data.t = obj.settings.grid;
            obj.data.T = length(obj.data.t);
            
            output = obj.data;
        end
        
        
        function optimize(obj)
            tic
            
            grid = (obj.data.t - obj.data.t(1)) / range(obj.data.t);        % scaled time grid
            grid_fine = linspace(0, 1, 201);                                % scaled optimized time grid
            grid_fs = (obj.settings.grid - obj.data.t(1)) / range(obj.data.t);
            obj.T_fs = length(grid_fs);
            
            obj.B = cell(1, obj.L);
            obj.dB = cell(1, obj.L);
            obj.ddB = cell(1, obj.L);
            obj.B_fine = cell(1, obj.L);
            obj.dB_fine = cell(1, obj.L);
            obj.data.basis_fs = cell(1, obj.L);
            obj.data.dbasis_fs = cell(1, obj.L);
            obj.penalty_ind = cell(1, obj.L);
            
            for k = 1:obj.L
                obj.B{k} = obj.bsplines{k}.basis(grid);                     % B-spline basis evaluated at t
                obj.dB{k} = obj.bsplines{k}.basis_diff(grid);               % corresponding first derivative
                obj.data.basis_fs{k} = obj.bsplines{k}.basis(grid_fs);      % save for later use
                obj.data.dbasis_fs{k} = obj.bsplines{k}.basis(grid_fs); 
                obj.B_fine{k} = obj.bsplines{k}.basis(grid_fine);           % save for smooth plotting
                obj.dB_fine{k} = obj.bsplines{k}.basis_diff(grid_fine);
                
                if obj.settings.penalization == "curvature"
                    obj.ddB{k} = obj.bsplines{k}.basis_ddiff(grid);         % corresponding second derivative
                    obj.penalty_ind{k} = obj.settings.penalized(1, k) <= obj.data.t ...
                                             & obj.data.t <= obj.settings.penalized(2, k);
                                         
                elseif obj.settings.penalization == "coefficients"
                    D = zeros(obj.bsplines{k}.card);
                    D = spdiags(repmat([-1 2 -1], obj.bsplines{k}.card, 1), [-1 0 1], D);
                    obj.ddB{k} = D(2:end-1, :)';
                    penalty_ind_t = obj.settings.penalized(1, k) <= obj.data.t ...
                                    & obj.data.t <= obj.settings.penalized(2, k);
                    penalty_ind_combined = movmin(any(obj.B{k} > 0 & penalty_ind_t, 2), 3);
                    obj.penalty_ind{k} = penalty_ind_combined(2:end-1)';
                end
            end

            obj.data.basis = obj.B;                                         % save for later use
            obj.data.dbasis = obj.dB;
            obj.data.ddbasis = obj.ddB;
            obj.data.penalty_ind = obj.penalty_ind;

            for iter = 1:obj.settings.niter                                 % FWLS spline coefficient estimates
                delta_old = obj.delta; 
                obj.update_coefficients();                     
                obj.update_variances();
                if eucl_rel(delta_old', obj.delta') < obj.settings.tol        % check convergence
                    break
                end
            end
            
            obj.data.lambda = obj.lambda;                                   % save penalty multiplier
            obj.data.variances_sm = obj.variances_sm;
            obj.data.variances = obj.variances_fs;                          % save variance estimates
            obj.fit();                                                      % compute splines
            
            obj.data.toc_sm = toc;
        end
        
                                            
        function update_coefficients(obj)                               % Spline coefficients from current variances
            if ~obj.settings.interactive && isempty(obj.settings.lambda)
                obj.lambda = obj.GCV_lambda();                              % minimize Generalized CV error
            end
            
            for i = 1:obj.N                                                 % W: correcting weights from variances
                for k = 1:obj.L                                             % penalty: lambda * L2(spline basis/coefficients)
                    W = diag(max(obj.variances_sm(:, k, i), 1e-7).^-1);     % delta_ik: LS estimate
                    penalty = obj.lambda(k) * obj.ddB{k} * diag(obj.penalty_ind{k}) * obj.ddB{k}';
                    obj.delta{k}(:, i) = svdinv(obj.B{k} * W * obj.B{k}' + penalty) ...
                                         * obj.B{k} * W * obj.data.traces(:, k, i);
                end
            end
        end
        
        
        function update_variances(obj)                                  % Estimate variances from current splines
            for k = 1:obj.L
                predicted = obj.B{k}' * reshape(obj.delta{k}, obj.bsplines{k}.card, obj.N);
                predicted_fs = obj.data.basis_fs{k}' * reshape(obj.delta{k}, obj.bsplines{k}.card, obj.N);
                design = [ones(obj.N * obj.T, 1) flatten(predicted).^2];
                design_fs = [ones(obj.N * obj.T_fs, 1) flatten(predicted_fs).^2];
                response = (flatten(predicted) - flatten(obj.data.traces(:, k, :))).^2;
                
                coefficients = Optimizer.optimize(design, response);        % find nonzero LS estimates for theta
                if sum(coefficients) == 0, coefficients(1) = mean(response); end
                                                                            % optimize further iteratively
                coefficients = Optimizer.newton_raphson(coefficients, predicted, obj.data.traces(:, k, :));
                if sum(coefficients) == 0, coefficients(1) = mean(response); end
                
                obj.sigma(k) = coefficients(1);                                % store opt and compute variances
                obj.tau(k) = coefficients(2);
                obj.variances_sm(:, k, :) = reshape(design * coefficients', obj.T, 1, obj.N);
                obj.variances_fs(:, k, :) = reshape(design_fs * coefficients', obj.T_fs, 1, obj.N);
            end
        end
        
        
        function lambda = GCV_lambda(obj)                               % Generalized CV minimizer
            options = optimoptions('fmincon', 'Display', 'off', 'StepTolerance', 1e-2);
            
            ncells = 3;
            lambda = zeros(obj.L, ncells);                                  % penalty multiplier
            for k = 1:obj.L
                for rep = 1:ncells
                    i = randi(obj.N);
                    lambda(k, rep) = exp(fmincon(@(log_lambda_k) obj.GCV(log_lambda_k, i, k), ...
                                                 log(obj.lambda(k)), [], [], [], [], -15, 15, [], options));
                end
            end
            
            lambda = mean(lambda, 2)';
        end
        
        
        function error = GCV(obj, log_lambda_k, i, k)
            Y = obj.data.traces(:, k, i);
            W = diag(max(obj.variances_sm(:, k, i), 1e-12).^-1);
            penalty = obj.ddB{k} * diag(obj.penalty_ind{k}) * obj.ddB{k}';
                                                                            % svdinv if necessary
            A = sqrt(W) * obj.B{k}' * ((obj.B{k} * W * obj.B{k}' + exp(log_lambda_k) * penalty) \ obj.B{k} * sqrt(W));
            error = obj.T * norm((eye(obj.T) - A) * Y)^2 / trace(eye(obj.T) - A)^2;
        end
        
        
        function fit(obj)                                               % Fitted splines from estimated coefficients
            obj.data.smoothed = zeros(obj.T_fs, obj.L, obj.N);              % smoothed data
            obj.data.dsmoothed = zeros(obj.T_fs, obj.L, obj.N);             % smoothed derivative
            obj.data.smoothed_fine = zeros(201, obj.L, obj.N);              % same for finer time grid
            obj.data.dsmoothed_fine = zeros(201, obj.L, obj.N);             % smoothed derivative
            
            for i = 1:obj.N                                                 % coefficients to spline interpolants
                for k = 1:obj.L
                    obj.data.smoothed(:, k, i) = obj.data.basis_fs{k}' * obj.delta{k}(:, i);
                    obj.data.dsmoothed(:, k, i) = obj.data.dbasis_fs{k}' * obj.delta{k}(:, i) / range(obj.data.t);
                    obj.data.smoothed_fine(:, k, i) = obj.B_fine{k}' * obj.delta{k}(:, i);
                    obj.data.dsmoothed_fine(:, k, i) = obj.dB_fine{k}' * obj.delta{k}(:, i) / range(obj.data.t);
                end
            end
            
            obj.data.smoothed = max(obj.data.smoothed, 1e-12);
        end
    end
end