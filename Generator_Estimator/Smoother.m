classdef Smoother < handle
    properties (Access = private)
        T                                                                   % number of time points
        T_fs                                                                % number of optimized time points
        L                                                                   % system dimension
        N                                                                   % number of cells
        penalty_ind
        
        
    end
    
    properties (Access = public)
        data                                                                % general data and output container
        system
        settings
        bspline                                                             % nested object controlling splines from B-splines
        
        lambda                                                              % penalty multiplier
        fine                                                                % fine grid for plotting
        X                                                                   % B-spline basis evaluated at t
        dX                                                                  % corresponding first derivative
        ddX                                                                 % corresponding second derivative
        X_fine
        dX_fine
        smoothed                                                            % smoothed data
        dsmoothed                                                           % smoothed derivative
        variances_sm
        variances_fs
        theta_fs
    end
    
    methods
        function obj = Smoother(data, system, settings)                 % Constructor
            obj.data = data;                                                % store experiment data
            obj.system = system;
            obj.settings = settings;
            obj.bspline = BSpline(settings.order, settings.knots);          % instantiate spline
            [obj.T, obj.L, obj.N] = size(obj.data.traces);                  % extract trajectory dimensions
            if ~isempty(settings.lambda)
                obj.lambda = reshape(settings.lambda, 1, []) .* ones(1, obj.L);              % initial penalty multiplier
            else
                obj.lambda = 1e-3 * ones(1, obj.L);
            end
            obj.variances_sm = zeros(obj.T, obj.L, obj.N);
            obj.variances_fs = zeros(length(settings.grid), obj.L, obj.N);  % measurement error variances
            obj.theta_fs = zeros(1, 2*obj.L);
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
            
            %%%
%             load('tempsmoothed.mat')
%             obj.data.smoothed(:, 2, :) = smoothed(:, 2, :);
%             obj.data.dsmoothed(:, 2, :) = dsmoothed(:, 2, :);
%             obj.data.smoothed_fine(:, 2, :) = smoothed_fine(:, 2, :);
%             obj.data.dsmoothed_fine(:, 2, :) = dsmoothed_fine(:, 2, :);
            %%%
            
            output = obj.data;
        end
        
        
        function optimize(obj)
            tic
            
            grid = (obj.data.t - obj.data.t(1)) / range(obj.data.t);        % scaled time grid
            obj.X = obj.bspline.basis(grid);                                % B-spline basis evaluated at t
            obj.dX = obj.bspline.basis_diff(grid);                          % corresponding first derivative

            obj.fine = linspace(0, 1, 201);
            obj.X_fine = obj.bspline.basis(obj.fine);                       % B-spline basis evaluated at t
            obj.dX_fine = obj.bspline.basis_diff(obj.fine);                 % corresponding first derivative

%             obj.ddX = obj.bspline.basis_ddiff(grid);                    % corresponding second derivative
%             obj.penalty_ind = obj.settings.linear(1) <= obj.data.t' & obj.data.t' <= obj.settings.linear(2);
%             disp(obj.penalty_ind)
            D = zeros(obj.bspline.card);
            D = spdiags(repmat([-1 2 -1], obj.bspline.card, 1), [-1 0 1], D);
            
            obj.ddX = D(2:end-1, :)';
%             obj.ddX = [D(:, 2) D(:, 2:end-1) D(:, end-1)];
            
            penalty_ind_t = obj.settings.linear(1) <= obj.data.t & obj.data.t <= obj.settings.linear(2);
            obj.penalty_ind = movmin(any(obj.X > 0 & penalty_ind_t, 2), 3);
            obj.penalty_ind = obj.penalty_ind(2:end-1);

            obj.data.basis = obj.X;                                         % save for later use
            obj.data.dbasis = obj.dX;
            obj.data.ddbasis = obj.ddX;
            obj.data.penalty_ind = obj.penalty_ind';
            
            grid = (obj.settings.grid - obj.data.t(1)) / range(obj.data.t); % scaled optimized time grid
            obj.T_fs = length(grid);
            obj.data.basis_fs = obj.bspline.basis(grid);                    % save for later use
            obj.data.dbasis_fs = obj.bspline.basis(grid); 

            variances = ones(obj.T, obj.L, obj.N);                          % time-dependent variances
            B_old = zeros(obj.bspline.card, obj.L, obj.N);                  % previous coefficients

            rel_tol = 1e-2;                                                 % relative error tolerance          
            for rep = 1:20                                                  % FWLS spline coefficient estimates
                B = obj.update_coefficients(variances);                     
                variances = obj.update_variances(B);
                if eucl_rel(B_old, B) < rel_tol, break, end                 % check convergence
                B_old = B;                                                  % update estimate
            end

            obj.data.lambda = obj.lambda;                                   % save penalty multiplier
            obj.data.variances_sm = variances;
            obj.data.variances = obj.variances_fs;                          % save variance estimates
            obj.fit(B);                                                     % compute splines
            
            obj.data.toc_sm = toc;
        end
        
                                            
        function B = update_coefficients(obj, variances)                % Spline coefficients from current variances
            if ~obj.settings.interactive && isempty(obj.settings.lambda)
                obj.lambda = obj.GCV_lambda(variances);                     % minimize Generalized CV error
%                 disp('New lambda values:');
%                 disp(obj.lambda)
            end
            
            B = zeros(obj.bspline.card, obj.L, obj.N);                      % B: spline coefficients
            for i = 1:obj.N                                                 % W: correcting weights from variances
                for k = 1:obj.L                                             % penalty: lambda * L2(spline basis)
                    W = diag(max(variances(:, k, i), 1e-7).^-1);            % B_ik: LS estimate
                    penalty = obj.lambda(k) * obj.ddX * diag(obj.penalty_ind) * obj.ddX';
                    B_ik = svdinv(obj.X * W * obj.X' + penalty) * obj.X * W * obj.data.traces(:, k, i);
                    B(:, k, i) = B_ik;
                end
            end
        end
        
        
        function variances = update_variances(obj, B)                   % Estimate variances from current splines
            predicted = reshape(obj.X' * reshape(B, obj.bspline.card, []), obj.T, obj.L, obj.N);
            R2 = (predicted - obj.data.traces).^2;
                                                                            % variant for optimized times
            predicted_fs = reshape(obj.data.basis_fs' * reshape(B, obj.bspline.card, []), obj.T_fs, obj.L, obj.N);

            theta = zeros(1, 2*obj.L);                                      % var_ik(t) = th1 + th2 * spline_ik(t)^2
            variances = zeros(obj.T, obj.L, obj.N);                         % time-dependent variances
            for k = 1:obj.L
                design = [ones(obj.N * obj.T, 1) flatten(predicted(:, k, :)).^2];
                response = flatten(R2(:, k, :));
                coefficients = Optimizer.optimize(design, response);        % find nonzero LS estimates for theta
                if sum(coefficients) == 0, coefficients(1) = mean(response); end
                                                                            % take two clipped Newton-Raphson steps
                coefficients = Optimizer.newton_raphson(coefficients, predicted(:, k, :), obj.data.traces(:, k, :));
                if sum(coefficients) == 0, coefficients(1) = mean(response); end
%                 coefficients(2) = coefficients(2) * 6;
                theta(2*(k-1)+[1 2]) = coefficients;                        % store opt and compute variances
                variances(:, k, :) = reshape(design * coefficients', obj.T, 1, obj.N);
                
                design_fs = [ones(obj.N * obj.T_fs, 1) flatten(predicted_fs(:, k, :)).^2];
                obj.variances_fs(:, k, :) = reshape(design_fs * coefficients', obj.T_fs, 1, obj.N);
            end
            
            obj.theta_fs = theta;
        end
        
        
        function lambda = GCV_lambda(obj, variances)                    % Generalized CV minimizer
            options = optimoptions('fmincon', 'Display', 'off', 'StepTolerance', 1e-2);
            
            ncells = 3;
            lambda = zeros(obj.L, ncells);                                  % penalty multiplier
            for k = 1:obj.L
                for rep = 1:ncells
                    i = randi(obj.N);
%                     Y = obj.data.traces(:, k, i);                           % corresponding trajectory
%                     W = diag(max(variances(:, k, i), 1e-12).^-1);            % correcting weights from variances
%       
%                     penalty = obj.ddX * diag(obj.penalty_ind) * obj.ddX';   % GCV
%                     A = @(la) sqrt(W) * obj.X' * svdinv(obj.X * W * obj.X' + la * penalty) ...
%                                                    * obj.X * sqrt(W);
%                     norm2 = @(x) sum(x.^2);
%                     auxiliary = @(A_la) obj.T * norm2((eye(obj.T) - A_la) * Y) / trace(eye(obj.T) - A_la)^2;
%                     objective = @(log_la) auxiliary(A(exp(log_la)));
% 
%                     lambda(k, rep) = exp(fmincon(objective, log(obj.lambda(k)), [], [], [], [], -15, 15, [], options));
                    
                    lambda(k, rep) = exp(fmincon(@(log_lambda_k) obj.GCV(log_lambda_k, i, k, variances), ...
                                                 log(obj.lambda(k)), [], [], [], [], -15, 15, [], options));
                end
            end
            
            lambda = mean(lambda, 2)';
        end
        
        
        function error = GCV(obj, log_lambda_k, i, k, variances)
            Y = obj.data.traces(:, k, i);
            W = diag(max(variances(:, k, i), 1e-12).^-1);
            penalty = obj.ddX * diag(obj.penalty_ind) * obj.ddX';
                                                                    % svdinv if necessary
            A = sqrt(W) * obj.X' * ((obj.X * W * obj.X' + exp(log_lambda_k) * penalty) \ obj.X * sqrt(W));
            error = obj.T * norm((eye(obj.T) - A) * Y)^2 / trace(eye(obj.T) - A)^2;
        end
        
        
%         function lambda = GCV_lambda(obj, variances)                    % Generalized CV minimizer
% %             options = optimoptions('fmincon', 'Display', 'off', 'OptimalityTolerance', 1e-3);            
%             options = optimoptions('fmincon', 'Display', 'off', 'StepTolerance', 1e-2);
%             
%             ncells = 3;
%             indices = sort(randperm(obj.N, ncells));
%             bound = 15 * ones(1, obj.L);
%             
%             lambda = exp(fmincon(@(log_lambda) obj.total_GCV(log_lambda, indices, variances), ...
%                                  log(obj.lambda), [], [], [], [], -bound, bound, [], options));
%         end
%         
%         
%         function error = total_GCV(obj, log_lambda, indices, variances)
%             error = 0;
%             for i = indices
%                 for k = 1:obj.L
%                     Y = obj.data.traces(:, k, i);
%                     W = diag(max(variances(:, k, i), 1e-12).^-1);
%                     penalty = obj.ddX * diag(obj.penalty_ind) * obj.ddX';
%                                                                             % svdinv if necessary
%                     A = sqrt(W) * obj.X' * ((obj.X * W * obj.X' + exp(log_lambda(k)) * penalty) \ obj.X * sqrt(W));
%                     error = error + obj.T * norm((eye(obj.T) - A) * Y)^2 / trace(eye(obj.T) - A)^2;
%                 end
%             end
%         end
        
        
        function fit(obj, B)                                            % Fitted splines from estimated coefficients         
            grid = (obj.settings.grid - obj.data.t(1)) / range(obj.data.t); % scaled optimized time grid
            obj.data.basis_fs = obj.bspline.basis(grid);                    % save for later use
            obj.data.dbasis_fs = obj.bspline.basis_diff(grid); 
            
            obj.data.smoothed = zeros(obj.T_fs, obj.L, obj.N);              % smoothed data
            obj.data.dsmoothed = zeros(obj.T_fs, obj.L, obj.N);             % smoothed derivative
            obj.data.smoothed_fine = zeros(201, obj.L, obj.N);              % same for finer time grid
            obj.data.dsmoothed_fine = zeros(201, obj.L, obj.N);             % smoothed derivative
            
            for i = 1:obj.N                                                 % coefficients to spline interpolants
                obj.data.smoothed(:, :, i) = obj.data.basis_fs' * B(:, :, i);
                obj.data.dsmoothed(:, :, i) = obj.data.dbasis_fs' * B(:, :, i) / range(obj.data.t);
                obj.data.smoothed_fine(:, :, i) = obj.X_fine' * B(:, :, i);
                obj.data.dsmoothed_fine(:, :, i) = obj.dX_fine' * B(:, :, i) / range(obj.data.t);
            end
            
            obj.data.smoothed = max(obj.data.smoothed, 1e-12);
        end
    end
end