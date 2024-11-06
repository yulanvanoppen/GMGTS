classdef SecondStage < handle
    properties (SetAccess = private)
        data                                                                % data aggregate struct
        system                                                              % nested object controlling the ODE system
        settings                                                            % hyperparameters and user input
        
        precision                                                           % cell parameter FS precision estimates
        beta                                                                % cell-specific parameter iterations
        b                                                                   % population parameter iterations
        D                                                                   % random effect covariance iterations
    end
    
    
    methods
        function obj = SecondStage(data, system, settings)              % Constructor
            obj.data = data;
            obj.system = system;
            obj.settings = settings;
            
            obj.beta = data.beta_fs;                                  % initial estimates from first stage results
%             obj.beta = data.beta_fs(obj.data.converged, :);
            obj.b = mean(data.beta_fs);                             % sample mean and covariance
            obj.D = cov(data.beta_fs);
        end
        
                                                                        % Main optimization function
        function output = optimize(obj)                                     % invert uncertainty estimates
            obj.estimate_precisions();                                      % correcting for near-zero variances
            
            for iter = 1:obj.settings.niter
                fixed = diag(obj.D(:, :, iter)) < max(diag(obj.D(:, :, iter))) / 1e6;
                Dinv_iter = zeros(size(obj.D(:, :, iter)));                   % invert current estimate of D
                Dinv_iter(~fixed, ~fixed) = svdinv(obj.D(~fixed, ~fixed, iter));
                Dinv_iter(fixed, fixed) = diag(Inf(1, sum(fixed)));
                                                                            % iterate until convergence:
                obj.E_step(iter, Dinv_iter)                                   % - refine cell-specific estimates
                obj.M_step(iter, Dinv_iter)                                   % - refine population estimates and covariances
                
                if eucl_rel(obj.b(iter, :), obj.b(iter+1, :)) < obj.settings.tol && ...
                   eucl_rel(obj.D(:, :, iter), obj.D(:, :, iter+1), true) < obj.settings.tol
                    break                                                   % break at convergence tolerance
                end
            end
            
            obj.extract_estimates();                                        % compute final predictions
            output = obj.data;                                              % return data struct extended with results
        end
        
        
        function E_step(obj, iter, Dinv_iter)                           % Expectation step of EM algorithm
            b_iter = obj.b(iter, :)';
            for i = 1:obj.data.N
                beta_fs_i = obj.beta(i, :, 1)';
                Cinv_i = obj.precision(:, :, i);
                
                fixed = isinf(diag(Cinv_i)) | isinf(diag(Dinv_iter));       % manually set for near-zero variability
                obj.beta(i, fixed, iter+1) = beta_fs_i(fixed);              % update estimate
                obj.beta(i, ~fixed, iter+1) = ((Cinv_i(~fixed, ~fixed) ...
                                                     + Dinv_iter(~fixed, ~fixed)) ...
                                                    \ (Cinv_i(~fixed, ~fixed) * beta_fs_i(~fixed) ...
                                                       + Dinv_iter(~fixed, ~fixed) * b_iter(~fixed)))';
            end
        end

        
        function M_step(obj, iter, Dinv_iter)                           % Maximization step of EM algorithm
            warning('off','MATLAB:singularMatrix')
            
            obj.b(iter+1, :) = mean(obj.beta(:, :, iter+1));                % update parameter mean
            b_iter = obj.b(iter+1, :)';

            summands = zeros(obj.system.P, obj.system.P, obj.data.N);       % collect terms of D
            for i = 1:obj.data.N
                beta_i = obj.beta(i, :, iter+1)';
                Cinv_i = obj.precision(:, :, i);
                
                fixed = isinf(diag(Cinv_i)) | isinf(diag(Dinv_iter));
                cov_term = zeros(size(Dinv_iter));                          % manually set for near-zero variability
                cov_term(~fixed, ~fixed) = svdinv(Cinv_i(~fixed, ~fixed) + Dinv_iter(~fixed, ~fixed));
                
                summands(:, :, i) = cov_term + (beta_i - b_iter) * (beta_i - b_iter)';
            end
            
            obj.D(:, :, iter+1) = mean(summands, 3);                        % update D
            warning('on')
        end

        
        function extract_estimates(obj)                                 % Collect final estimates and predictions
            obj.data.precision = obj.precision;
            obj.data.beta_ss = max(0, obj.beta(:, :, end));                % cell parametes
            obj.data.b_est = max(0, obj.b(end, :));                         % population means
            obj.data.D_est = obj.D(:, :, end);                              % random effect covariance
                                                                            % fitted trajectories for empirical Bayes estimators
            obj.data.fitted2 = obj.system.integrate(obj.data.beta_ss, obj.data, obj.data.t_fine);
            
                                                                            % parameter population mean trajectory
            obj.data.population = obj.system.integrate(obj.data.b_est, obj.data, obj.data.t_fine);
        end
        
        
        function estimate_precisions(obj)                               % Inverse covariance matrix estimates
            for i = 1:obj.data.N
                covariance = obj.data.varbeta(:, :, i);
                
                fixed = diag(covariance) < max(diag(covariance)) / 1e6;     % manual correction for near-zero uncertainty
                obj.precision(~fixed, ~fixed, i) = svdinv(covariance(~fixed, ~fixed));
                obj.precision(fixed, fixed, i) = diag(Inf(1, sum(fixed)));
            end
        end
    end
end