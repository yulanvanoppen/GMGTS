% START FROM TRUE CELL PARAMETERS
% LOG-LIKELIHOOD VALUES
%


classdef FirstStageGTS < handle
    properties (SetAccess = private)
        data                                                                % data aggregate struct
        system                                                              % nested object controlling the ODE system
        settings                                                            % weights and loop controls
        
        T                                                                   % number of time points
        L                                                                   % number of observed states
        N                                                                   % population size
        
        beta_fs                                                             % cell-specific parameter estimates 
        varbeta                                                             % uncertainty estimates of beta
        sigma2                                                              % additive measurement error variance 
        tau2                                                                % multiplicative measurement error variance 
        variances_fs                                                        % measurement error variances along time grid
        
        convergence_steps                                                   % relative iteration steps
        not_converged                                                       % indices of cells that have not yet converged
        
        fitted_fs                                                           % predicted states and gradients
        dfitted_fs
    end
    
    
    methods                                                 
        function obj = FirstStageGTS(data, system, settings)            % Constructor
            obj.data = data;                                                % experimental and smoothing data
            obj.system = system;                                            % ODE system
            obj.settings = settings;                                        % hyperparameters
            
            [obj.T, obj.L, obj.N] = size(obj.data.traces);                  % extract from trajectory dimensions
            
            obj.varbeta = repmat(eye(system.P), 1, 1, obj.N);
            obj.sigma2 = zeros(1, obj.L);
            obj.tau2 = zeros(1, obj.L);
            obj.variances_fs = ones(obj.T, obj.L, obj.N);
            
            obj.convergence_steps = ones(1, obj.N);                         % ensure no cells are considered converged
            obj.not_converged = 1:obj.N;
            
            obj.initialize();                                               % numerically optimize population average
        end
        
        
        function output = optimize(obj)                                 % Main optimization function
            disp('Classical GTS: Trajectory matching')
            
            for iter = 1:obj.settings.niter
                beta_old = obj.beta_fs;
                
                obj.update_parameters(iter);                                % update (cell-specific) parameter estimates
                obj.estimate_variances();                                   % update variance estimates
                
                if obj.system.P > 1                                         % compute relative iteration steps
                    obj.convergence_steps = vecnorm((beta_old - obj.beta_fs)') ./ vecnorm(beta_old');
                else
                    obj.convergence_steps = abs(beta_old' - obj.beta_fs') ./ abs(beta_old');
                end                                                         % cells that have not converged 
                obj.not_converged = find(obj.convergence_steps >= obj.settings.tol);
                
                fprintf('%3d (%.1e, %.1e, %.1e)\n', iter, ...               % print iteration step quantiles .5, .9, 1.0
                        median(obj.convergence_steps), quantile(obj.convergence_steps, .9), max(obj.convergence_steps));
                
                if quantile(obj.convergence_steps, .9) < obj.settings.tol, break, end
            end
            
            obj.uncertainties_beta();                                       % compute estimate uncertainties
            
            obj.save_estimates();                                           % save estimates and fitted trajectories
            output = obj.data;                                              % return data appended with FS results
        end
        
        
        function initialize(obj)                                        % Initial numerical optimization on population average
            beta_init = Optimization.least_squares(@obj.squares_sum, obj.settings.lb, obj.settings.ub, obj.settings.nstart);
            obj.beta_fs = repmat(beta_init, obj.N, 1);                      % copy to each cell
        end
        
        
        function update_parameters(obj, iter)                           % Update (cell-specific) parameter estimates
            for i = 1:obj.N
                SS = @(beta) obj.squares_sum(beta, obj.data.traces(obj.settings.optindices, :, i), ...
                                                   obj.variances_fs(obj.settings.optindices, :, i));
                                               
                if iter == 1                                                % multistart every cell (for TCS)
                    obj.beta_fs(i, :) = Optimization.least_squares(SS, obj.settings.lb, obj.settings.ub, 5, []);
                else
                    obj.beta_fs(i, :) = Optimization.least_squares(SS, obj.settings.lb, obj.settings.ub, 1, obj.beta_fs(i, :));
                end
            end
            obj.fitted_fs = obj.system.integrate(obj.beta_fs, obj.data);    % compute predicted trajectories
            if obj.settings.positive
                obj.fitted_fs = max(1e-12, obj.fitted_fs);                  % force positive
            end
        end
        
        
        function ss = squares_sum(obj, beta, y, variances)              % Weighted sum of squared differences
            ss = Inf;
            try                                                 
                solution = obj.system.integrate(beta, obj.data);            % compute fitted trajectories
                solution = solution(obj.settings.optindices, obj.data.observed, :);
                
                if nargin == 2                                              % sum of squares on observed states
                    ss = sum((solution - obj.data.traces(obj.settings.optindices, :, :)).^2 ...
                             ./ mean(obj.data.traces(obj.settings.optindices, :, :), [1 3]).^2, 'all');
                                                                            % add prior term (if any)
                    ss = ss + (beta' - obj.settings.prior.mean)' * (obj.settings.prior.prec * obj.settings.prior.mult) ...
                                                                 * (beta' - obj.settings.prior.mean);
                else
                    ss = sum((solution - y).^2 ./ variances, 'all');        % sum of squares on observed states
                                                                            % add prior term (if any)
                    ss = ss + (beta' - obj.settings.prior.mean)' * obj.settings.prior.prec ...
                                                                 * (beta' - obj.settings.prior.mean);
                end
            catch ME
                disp(ME)
            end
        end
        
        
        function estimate_variances(obj)                                % Update measurement error variances
            for k = 1:obj.L
                predicted = obj.fitted_fs(:, obj.data.observed(k), :);      % fitted trajectories
                design = [ones(obj.N*obj.T, 1) flatten(predicted).^2];      % columns for additive and multiplicative noise
                response = flatten(predicted - obj.data.traces(:, k, :)).^2;
                
                coefficients = lsqnonneg(design, response)';                % initialize nonzero LS estimates for noise parameters
                if sum(coefficients) == 0, coefficients(1) = mean(response); end
                                                                            % optimize further iteratively
                coefficients = Optimization.noise_parameters(coefficients, predicted, obj.data.traces(:, k, :));
                if sum(coefficients) == 0, coefficients(1) = mean(response); end
                
                obj.sigma2(k) = coefficients(1);                            % store optimum and compute variances
                obj.tau2(k) = coefficients(2);
                
                obj.variances_fs(:, k, :) = reshape(design * coefficients', obj.T, 1, obj.N);
            end
        end
        
        
        function save_estimates(obj)                                    % Extract results 
            obj.data.beta_fs = obj.beta_fs;
            obj.data.fitted_fs = obj.fitted_fs;           
            obj.data.varbeta = obj.varbeta;
            obj.data.variances_fs = obj.variances_fs;
            obj.data.convergence_steps = obj.convergence_steps;
            obj.data.converged = setdiff(1:obj.N, obj.not_converged);       % integrate along finer time grid
            obj.data.fitted_fs_fine = obj.system.integrate(obj.beta_fs, obj.data, obj.data.t_fine);
            obj.data.dfitted_fs_fine = obj.system.rhs(obj.data.fitted_fs_fine, obj.data.t_fine, obj.beta_fs);
            if obj.settings.positive
                obj.data.fitted_fs_fine = max(1e-12, obj.data.fitted_fs_fine);
            end
            if obj.settings.lognormal, obj.lognormal_approximation(), end
        end
        
        
        function lognormal_approximation(obj)
            obj.data.beta_lnorm = zeros(size(obj.beta_fs));
            obj.data.varbeta_lnorm = zeros(size(obj.varbeta));
            for i = 1:obj.N
                beta_i = obj.beta_fs(i, :)';
                varbeta_i = obj.varbeta(:, :, i);
                                                                            % moment matching
%                 obj.data.varbeta_lnorm(:, :, i) = log(1 + varbeta_i ./ (beta_i*beta_i'));
%                 obj.data.beta_lnorm(i, :) = log(beta_i') - .5 * diag(obj.data.varbeta_lnorm(:, :, i))';
                
                                                                            % unbiased log
                diag_i = diag(log(.5 + sqrt(.25 + diag(diag(varbeta_i)./beta_i.^2))));
                full_i = log(varbeta_i ./ (beta_i * beta_i') .* exp(-.5 * (diag_i + diag_i')) + 1);
                obj.data.varbeta_lnorm(:, :, i) = full_i - diag(diag_i - diag(full_i));
                obj.data.beta_lnorm(i, :) = log(beta_i');
                if ~isreal(obj.data.beta_lnorm(i, :)) || ~isreal(obj.data.varbeta_lnorm(:, :, i))
                    obj.data.convergence_steps(i) = 1;
                end
            end
        end
        
        
        function uncertainties_beta(obj)                                % Parameter uncertainties
            beta_permuted = permute(obj.beta_fs, [3 2 1]);                  % vectorized finite difference approximations
            epsilon = max(1e-8, beta_permuted * .001);                      % of solution and RHS parameter sensitivities
            beta_pm_eps = beta_permuted + epsilon .* kron([1; -1], eye(obj.system.P));

            traces_pm_eps = zeros(obj.T, obj.system.K, obj.system.P, 2, obj.N);
            for i = 1:obj.N
                for p = 1:obj.system.P                                      % numerically integrate at parameter component pm eps
                    traces_pm_eps(:, :, p, 1, i) = obj.system.integrate(beta_pm_eps(p, :, i), obj.data, obj.data.t, 1e-4);
                    traces_pm_eps(:, :, p, 2, i) = obj.system.integrate(beta_pm_eps(p+end/2, :, i), obj.data, obj.data.t, 1e-4);
                end
            end

            eps_denom = 1./reshape(2*epsilon, 1, 1, obj.system.P, 1, obj.N);% finite difference approximations
            dF_dbeta_all = permute((traces_pm_eps(:, obj.data.observed, :, 1, :) ...
                                    - traces_pm_eps(:, obj.data.observed, :, 2, :)) .* eps_denom, [1:3 5 4]);
                                
            for i = 1:obj.N
                gradient = permute(dF_dbeta_all(:, :, :, i), [3 4 1 2]);    % scale by variances
                scaled_gradient = gradient ./ reshape(sqrt(obj.variances_fs(:, :, i)), ...
                                                      1, 1, obj.T, obj.L);  % inner products
                gradient2 = pagemtimes(scaled_gradient, 'none', scaled_gradient, 'transpose');
                precision = sum(gradient2(:, :, :, :), [3 4]);              % sum across time points and states
                obj.varbeta(:, :, i) = tryinv(precision);
            end
        end
    end
end












