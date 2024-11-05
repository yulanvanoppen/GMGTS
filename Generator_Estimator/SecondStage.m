classdef SecondStage < handle
    properties (SetAccess = private)
        data                                                                % data aggregate struct
        system                                                              % nested object controlling the ODE system
        settings                                                            % hyperparameters and user input
        
        precision                                                           % cell parameter FS precision estimates
        beta_cells                                                          % cell-specific parameter iterations
        beta_mean                                                           % population parameter iterations
        D                                                                   % random effect covariance iterations
    end
    
    
    methods
        function obj = SecondStage(data, system, settings)              % Constructor
            obj.data = data;
            obj.system = system;
            obj.settings = settings;
            
            obj.beta_cells = data.beta_fs;                                  % initial estimates from first stage results
%             obj.beta_cells = data.beta_fs(obj.data.converged, :);
            obj.beta_mean = mean(data.beta_fs);                             % sample mean and covariance
            obj.D = cov(data.beta_fs);
        end
        
                                                                        % Main optimization function
        function output = optimize(obj)                                     % invert uncertainty estimates
            obj.estimate_precisions();                                      % correcting for near-zero variances
            
            for rep = 1:obj.settings.nrep
                fixed = diag(obj.D(:, :, rep)) < max(diag(obj.D(:, :, rep))) / 1e6;
                Dinv_rep = zeros(size(obj.D(:, :, rep)));                   % invert current estimate of D
                Dinv_rep(~fixed, ~fixed) = svdinv(obj.D(~fixed, ~fixed, rep));
                Dinv_rep(fixed, fixed) = diag(Inf(1, sum(fixed)));
                                                                            % iterate until convergence:
                obj.E_step(rep, Dinv_rep)                                   % - refine cell-specific estimates
                obj.M_step(rep, Dinv_rep)                                   % - refine population estimates and covariances
                
                if eucl_rel(obj.beta_mean(rep, :), obj.beta_mean(rep+1, :)) < obj.settings.tol && ...
                   eucl_rel(obj.D(:, :, rep), obj.D(:, :, rep+1), true) < obj.settings.tol
                    break                                                   % break at convergence tolerance
                end
            end
            
            obj.extract_estimates();                                        % compute final predictions
            output = obj.data;                                              % return data struct extended with results
        end
        
                                            % Expectation step of EM algorithm
        function E_step(obj, rep, Dinv_rep)     % abbreviate for clarity
            beta_rep = obj.beta_mean(rep, :)';

                                                % copy to retain non-varying parameters
            obj.beta_cells(:, :, rep+1) = obj.beta_cells(:, :, rep);

            for i = 1:obj.data.N                     % cell i's initial beta and Ci^-1
                init_i = obj.beta_cells(i, :, 1)';
                Cinv_i = obj.precision(:, :, i);
                
                fixed = isinf(diag(Cinv_i)) | isinf(diag(Dinv_rep));
                
                                                % update estimate
                obj.beta_cells(i, fixed, rep+1) = init_i(fixed);
                obj.beta_cells(i, ~fixed, rep+1) = ((Cinv_i(~fixed, ~fixed) ...
                                                     + Dinv_rep(~fixed, ~fixed)) ...
                                                    \ (Cinv_i(~fixed, ~fixed) * init_i(~fixed) ...
                                                       + Dinv_rep(~fixed, ~fixed) * beta_rep(~fixed)))';
            end
        end

                                            % Maximization step of EM algorithm
        function M_step(obj, rep, Dinv_rep)
            warning('off','MATLAB:singularMatrix')
                                                % update population parameters and abbreviate
            obj.beta_mean(rep+1, :) = obj.beta_mean(rep, :);
            obj.beta_mean(rep+1, :) = mean(obj.beta_cells(:, :, rep+1));
            beta_rep = obj.beta_mean(rep+1, :)';

                                                % collect terms of D
            summands = zeros(obj.system.P, obj.system.P, obj.data.N);

            for i = 1:obj.data.N                     % abbreviate cell i's (current) beta and C^-1
                beta_i = obj.beta_cells(i, :, rep+1)';
                Cinv_i = obj.precision(:, :, i);
                
                certain = isinf(diag(Cinv_i)) | isinf(diag(Dinv_rep));
                cov_term = zeros(size(Dinv_rep));
                cov_term(~certain, ~certain) = svdinv(Cinv_i(~certain, ~certain) + Dinv_rep(~certain, ~certain));
                
                summands(:, :, i) = cov_term + (beta_i - beta_rep) * (beta_i - beta_rep)';
            end
                                                % update D
            obj.D(:, :, rep+1) = mean(summands, 3);
            warning('on')
        end

                                            % Extract results
        function extract_estimates(obj)         % collect final estimates:
            obj.data.precision = obj.precision;
            
            obj.data.beta_est = max(0, obj.beta_cells(:, :, end));  % - cell parametes
            obj.data.b_est = max(0, obj.beta_mean(end, :));         % - population means


            obj.data.D_GTS = obj.D(:, :, end);                      % - random effect covariance
                                               % compute final fitted cell trajectories
            obj.data.fitted2 = obj.system.integrate(obj.data.beta_est, obj.data, obj.data.t_fine);
            
                                               % compute population mean trajectory
            obj.data.population = obj.system.integrate(obj.data.b_est, obj.data, obj.data.t_fine);
        end
        
        
        function estimate_precisions(obj)     % inverse covariance matrix estimates
            for i = 1:obj.data.N
                covariance = obj.data.variances_beta_fs(:, :, i);
                
                fixed = diag(covariance) < max(diag(covariance)) / 1e6;
                obj.precision(~fixed, ~fixed, i) = svdinv(covariance(~fixed, ~fixed));
                obj.precision(fixed, fixed, i) = diag(Inf(1, sum(fixed)));
            end
        end
    end
end