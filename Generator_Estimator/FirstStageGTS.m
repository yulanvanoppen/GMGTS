classdef FirstStageGTS < handle
    properties (Access = private)
        data                                                                % data aggregate struct
        system                                                              % nested object controlling the ODE system
        settings                                                            % weights and loop controls
        
        T                                                                   % number of time points
        L                                                                   % number of observed states
        N                                                                   % population size
    end
    
    properties (Access = public)
        beta_fs                                                             % cell-specific parameter estimates
        sigma_fs                                                            % estimated measurement error magnitude   
        variances_fs
        not_converged
        theta_fs
        
        fitted_fs                                                           % predictions + derivatives
        dfitted_fs
    end
    
    methods                                                             % Constructor
        function obj = FirstStageGTS(data, system, settings)
            obj.data = data;
            obj.system = system;
            obj.settings = settings;
            
            [obj.T, obj.L, obj.N] = size(obj.data.traces);   % extract from trajectory dimensions
            obj.beta_fs = repmat(settings.initial, obj.N, 1);               % repeat initial estimate for each cell
            obj.sigma_fs = ones(2, obj.system.K);
            obj.variances_fs = ones(obj.T, obj.L, obj.N);
            obj.theta_fs = zeros(1, 2*obj.L);
            
            obj.initialize();
        end
        
        
        function output = optimize(obj)                                 % Main optimization function
            disp('Classical GTS: Trajectory matching')
            
            for rep = 1:obj.settings.nrep
                beta_old = obj.beta_fs;
                
                obj.update_parameters();                                    % update (cell-specific) parameter estimates
                obj.update_covariances();                                   % update variance estimates
                
                if obj.system.P > 1
                    solnorms = vecnorm((beta_old - obj.beta_fs)') ./ vecnorm(beta_old');
                else
                    solnorms = abs(beta_old' - obj.beta_fs') ./ abs(beta_old');
                end
                obj.not_converged = find(solnorms >= obj.settings.tol);
                
                fprintf('%3d (%.1e, %.1e, %.1e)\n', rep, ...
                        median(solnorms), quantile(solnorms, .9), max(solnorms));
                
                if quantile(solnorms, .9) < obj.settings.tol, break, end
            end
            
            obj.beta_fs(obj.beta_fs < 1e-15) = 0;                           % polish final estimates
            obj.data.beta_fs = obj.beta_fs;                                 % store results
            obj.data.fitted_fs = obj.fitted_fs;
            obj.data.sigma_fs = obj.sigma_fs;
            obj.data.variances_fs = obj.variances_fs;
            output = obj.data;                                              % return expanded data container
        end
        
        
        function initialize(obj)
            options = optimoptions('fmincon', 'Display', 'off', 'OptimalityTolerance', obj.settings.tol);
            value = Inf;
                                                                            % multistart on first run
            for start = 1:obj.settings.nstart
                logb0 = rand(1, obj.system.P) .* (log(obj.settings.ub) - log(obj.settings.lb)) + log(obj.settings.lb);
                [opt_new, value_new] = fmincon(@(logb0) obj.squares_sum_all(exp(logb0)), logb0, [], [], [], [], ...
                                               log(obj.settings.lb), log(obj.settings.ub), [], options);
                if value_new < value, value = value_new; opt = opt_new; end
            end
            disp("Initial beta:")
            disp(exp(opt))
            obj.beta_fs = repmat(exp(opt), obj.N, 1) .* (1 + obj.settings.perturbation * randn(obj.N, obj.system.P));
%             obj.beta_fs = obj.data.beta;
        end
        
        
        function ss = squares_sum_all(obj, b)              % Weighted um of squared differences
            ss = Inf;
            try
                solution = obj.system.integrate(b, obj.data);
                solution = solution(obj.settings.optindices, obj.data.observed, :);
                ss = sum((solution - obj.data.traces(obj.settings.optindices, :, :)).^2 ...
                         ./ mean(obj.data.traces(obj.settings.optindices, :, :), [1 3]).^2, 'all');
                ss = ss + (b' - obj.settings.prior.mean)' * (obj.settings.prior.prec * obj.settings.prior.mult) ...
                                                          * (b' - obj.settings.prior.mean);
            catch ME
                disp(ME)
            end
        end
        
        
        function update_parameters(obj)                                                  % population estimate
            estimates = obj.beta_fs;
            for i = 1:obj.N                                                 % update estimates
                estimates(i, :) = obj.minimize_SS(estimates(i, :), i);
            end
            obj.beta_fs = estimates;
            obj.fitted_fs = obj.system.integrate(obj.beta_fs, obj.data);    % compute predicted trajectories
        end
        
        
        function update_covariances(obj)                                % Update measurement error variances
            predicted = obj.fitted_fs(:, obj.data.observed, :);
            R2 = (predicted - obj.data.traces).^2;                          % squared differences

            theta = zeros(1, 2*obj.L);                                      % var_ik(t) = th1 + th2 * spline_ik(t)^2
            obj.variances_fs = ones(obj.T, obj.L, obj.N);                   % time-dependent variances
            for k = 1:obj.L
                design = [ones(obj.N * obj.T, 1) flatten(predicted(:, k, :)).^2];
                response = flatten(R2(:, k, :));
                coefficients = Optimizer.optimize(design, response);        % find nonzero LS estimates for theta
                if sum(coefficients) == 0, coefficients(1) = mean(response); end
                coefficients = coefficients + 1e-7;                         % take two clipped Newton-Raphson steps
                coefficients = Optimizer.newton_raphson(coefficients, predicted(:, k, :), obj.data.traces(:, k, :));
                theta(2*(k-1)+[1 2]) = coefficients;                        % store opt and compute variances
                obj.variances_fs(:, k, :) = reshape(design * coefficients', obj.T, 1, obj.N);
            end
            
            obj.theta_fs = theta;
        end
        
                                                                        % Minimize weighted sum of squared differences
        function opt = minimize_SS(obj, init, index)
            options = optimoptions('fmincon', 'Display', 'off', 'OptimalityTolerance', 1e-6);
            SS = @(beta) obj.squares_sum(beta, obj.data.traces(obj.settings.optindices, :, index), ...
                                                               obj.variances_fs(obj.settings.optindices, :, index));
            opt = zeros(size(init));
            value = Inf;                                                    % multistart on first run
            for start = 1 %:(1 + (obj.settings.nstart-1) * (all(~diff(obj.beta_fs), 'all')))
                b0 = init;
%                 if all(~diff(obj.beta_fs), 'all')
%                     b0 = exp(rand(1, obj.system.P) .* (log(obj.settings.ub) - log(obj.settings.lb)) + log(obj.settings.lb));
%                 end
                [opt_new, value_new] = fmincon(SS, b0, [], [], [], [], obj.settings.lb, obj.settings.ub, [], options);
                if value_new < value, value = value_new; opt = opt_new; end
            end
        end
        
        
        function ss = squares_sum(obj, beta, y, variances)              % Weighted um of squared differences
            ss = Inf;
            try
                solution = obj.system.integrate(beta, obj.data);
                solution = solution(obj.settings.optindices, obj.data.observed, :);
                ss = sum((solution - y).^2 ./ variances, 'all');
                ss = ss + (beta' - obj.settings.prior.mean)' * obj.settings.prior.prec * (beta' - obj.settings.prior.mean);
            catch ME
                disp(ME)
            end
        end
    end
end












