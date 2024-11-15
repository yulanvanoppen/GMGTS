classdef SecondStageGTS < handle
    properties (Access = private)
        data
        system                                                              % nested object controlling the ODE system
        settings
        
        varying                                                             % indices of parameters that may vary among cells
        nrep                                                                % maximum number of iterations
        tol                                                                 % desired error tolerance
        
        T                                                                   % number of time points
        L                                                                   % number of observed states
        N                                                                   % population size
        
        precision                                                           % cell parameter FS precision estimates
        beta                                                                % cell-specific parameter iterations
        b                                                                   % population parameter iterations
        D                                                                   % random effect covariance iterations
    end
    
    properties (SetAccess = public)
        beta_est                                                            % final cell-specific parameter estimates
        b_est                                                               % final population parameter estimates
        
        D_NTS                                                               % naive random effect covariance estimate
        D_GTS                                                               % final random effect covariance estimate
        
        fitted2                                                             % fitted cell trajectories
        population                                                          % fitted population trajectories
    end
    
    methods
        function obj = SecondStage(data, system, settings)              % Constructor
            obj.data = data;
            obj.system = system;
            obj.settings = settings;
            
            obj.beta = data.beta_fs;                                        % initial estimates from first stage results
            if settings.lognormal
                obj.beta = data.beta_lnorm;
            end                           
%             obj.beta = data.beta_fs(obj.data.converged, :);
            obj.b = mean(obj.beta);                                         % sample mean and covariance
            obj.D = cov(obj.beta);
        end
        

                                                                        % Main optimization function
        function output = optimize(obj)                                     % invert uncertainty estimates
            obj.estimate_precisions();                                      % correcting for near-zero variances
            
            for iter = 1:obj.settings.niter
                fixed = diag(obj.D(:, :, iter)) < max(diag(obj.D(:, :, iter))) / 1e6;
                Dinv_iter = zeros(size(obj.D(:, :, iter)));                 % invert current estimate of D
                Dinv_iter(~fixed, ~fixed) = svdinv(obj.D(~fixed, ~fixed, iter));
                Dinv_iter(fixed, fixed) = diag(Inf(1, sum(fixed)));
                                                                            % iterate until convergence:
                obj.E_step(iter, Dinv_iter)                                 % - refine cell-specific estimates
                obj.M_step(iter, Dinv_iter)                                 % - refine population estimates and covariances
                
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
            if obj.settings.lognormal
                obj.data.beta_ss = exp(obj.beta(:, :, end));
                obj.data.b_est = obj.b(end, :);                             % parameter population mean trajectory
                obj.data.population = obj.system.integrate(exp(obj.data.b_est), obj.data, obj.data.t_fine);
            else
                obj.data.beta_ss = max(0, obj.beta(:, :, end));             % cell parametes
                obj.data.b_est = max(0, obj.b(end, :));                     % population means
                obj.data.population = obj.system.integrate(obj.data.b_est, obj.data, obj.data.t_fine);
            end
            obj.data.D_est = obj.D(:, :, end);                              % random effect covariance
            obj.data.lognormal = obj.settings.lognormal;
                                                                            % fitted trajectories for empirical Bayes estimators
            obj.data.fitted2 = obj.system.integrate(obj.data.beta_ss, obj.data, obj.data.t_fine);
        end
        
        
        function estimate_precisions(obj)                               % Fischer information estimates
            obj.precision = zeros(length(obj.varying), length(obj.varying), obj.N);
           
            for i = 1:obj.N
                beta_i = obj.beta_cells(i, :);
                gradient = zeros(length(obj.varying), 1, obj.T, obj.L);
                for p = 1:length(obj.varying)
                    epsilon = max(1e-8, beta_i(obj.varying(p)) * .01);      % approximate using finite differences
                    beta_plus_eps = beta_i + epsilon * (obj.varying(p)==(1:obj.P));
                    beta_minus_eps = beta_i - epsilon * (obj.varying(p)==(1:obj.P));
                    traces_plus_eps = obj.system.integrate(beta_plus_eps, obj.data);
                    traces_minus_eps = obj.system.integrate(beta_minus_eps, obj.data);
                    full_gradient = (traces_plus_eps - traces_minus_eps) / (2*epsilon);
                    gradient(p, 1, :, :) = full_gradient(:, obj.data.observed, :);
                end                                                         % scale by variances
                scaled_gradient = gradient ./ reshape(sqrt(obj.data.variances_fs(:, :, i)), ...
                                                      1, 1, obj.T, obj.L);  % inner products
                gradient2 = pagemtimes(scaled_gradient, 'none', scaled_gradient, 'transpose');
                obj.precision(:, :, i) = sum(gradient2(:, :, :, :), [3 4]); % sum across time points and states
            end
        end
        
    end
end