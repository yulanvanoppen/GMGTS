classdef SecondStageGTS < handle
    properties (Access = private)
        data
        system                                                              % nested object controlling the ODE system
        settings
        
        varying                                                             % indices of parameters that may vary among cells
        nrep                                                                % maximum number of iterations
        tol                                                                 % desired error tolerance
        
        T                                                                   % number of time points
        K                                                                   % system dimension
        L                                                                   % number of observed states
        N                                                                   % population size
        P                                                                   % parameter vector length
        
        precision                                                           % asymptotic precision estimate
        beta_cells                                                          % cell-specific parameter iterations
        beta_mean                                                           % population parameter iterations
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
        function obj = SecondStageGTS(data, system, settings)           % Constructor
            obj.data = data;
            obj.system = system;
            obj.settings = settings;
            
            obj.varying = data.varying;
            
            [obj.T, obj.L, obj.N] = size(obj.data.traces);
            obj.K = system.K;
            obj.P = system.P;
            
            obj.beta_cells = data.beta_fs;
            obj.beta_mean = mean(data.beta_fs);
            obj.D = cov(data.beta_fs(:, obj.varying));
        end
        

        function output = optimize(obj)                                 % Main optimization function
            obj.estimate_precisions();                                      % approximate first-stage precisions
            
            for rep = 1:obj.settings.nrep
                certain = diag(obj.D(:, :, rep)) < 1e-7;                    % inversion of possibly singular matrix
                Dinv_rep = zeros(size(obj.D(:, :, rep)));
                Dinv_rep(~certain, ~certain) = inv(obj.D(~certain, ~certain, rep));
                Dinv_rep(certain, certain) = diag(Inf(1, sum(certain)));
                
                                                                            % iterate until convergence:
                obj.E_step(rep, Dinv_rep)                                   %   - refine cell-specific estimates
                obj.M_step(rep, Dinv_rep)                                   %   - refine population estimates and covariances
                
                                                                            % break upon convergence
                if eucl_rel(obj.beta_mean(rep, :), obj.beta_mean(rep+1, :)) < obj.settings.tol && ...
                   eucl_rel(obj.D(:, :, rep), obj.D(:, :, rep+1), true) < obj.settings.tol
                    break
                end
            end
            
            fprintf('\n')
            obj.extract_estimates();                                        % collect final estimates in obj.data
            output = obj.data;                                              % and return
        end
        

        function E_step(obj, rep, Dinv_rep)                             % Expectation step of EM algorithm
            beta_rep = obj.beta_mean(rep, obj.varying)';                    % abbreviate for clarity
            obj.beta_cells(:, :, rep+1) = obj.beta_cells(:, :, rep);        % copy to retain non-varying parameters

            for i = 1:obj.N
                init_i = obj.beta_cells(i, obj.varying, 1)';                % cell i's initial beta and Ci^-1
                Cinv_i = obj.precision(:, :, i);
                
                certain = isinf(diag(Cinv_i)) | isinf(diag(Dinv_rep));      % update estimate
                obj.beta_cells(i, obj.varying(certain), rep+1) = init_i(certain);
                obj.beta_cells(i, obj.varying(~certain), rep+1) = ((Cinv_i(~certain, ~certain) ...
                                                                    + Dinv_rep(~certain, ~certain)) ...
                                                                 \ (Cinv_i(~certain, ~certain) * init_i(~certain) ...
                                                                    + Dinv_rep(~certain, ~certain) * beta_rep(~certain)))';
            end
        end


        function M_step(obj, rep, Dinv_rep)                             % Maximization step of EM algorithm
            warning('off','MATLAB:singularMatrix')

            obj.beta_mean(rep+1, :) = obj.beta_mean(rep, :);                % update population parameters and abbreviate
            obj.beta_mean(rep+1, :) = mean(obj.beta_cells(:, :, rep+1));
            beta_rep = obj.beta_mean(rep+1, obj.varying)';
                                                                            % collect terms of D
            summands = zeros(length(obj.varying), length(obj.varying), obj.N);
            for i = 1:obj.N
                beta_i = obj.beta_cells(i, obj.varying, rep+1)';            % abbreviate cell i's (current) beta and C^-1
                Cinv_i = obj.precision(:, :, i);
                
                certain = isinf(diag(Cinv_i)) | isinf(diag(Dinv_rep));      % inversion of possibly singular matrix
                cov_term = zeros(size(Dinv_rep));
                cov_term(~certain, ~certain) = inv(Cinv_i(~certain, ~certain) + Dinv_rep(~certain, ~certain));
                
                summands(:, :, i) = cov_term + (beta_i - beta_rep) * (beta_i - beta_rep)';
            end

            obj.D(:, :, rep+1) = mean(summands, 3);                         % update D
            warning('on')
        end


        function extract_estimates(obj)                                 % Extract results
            obj.data.precision = obj.precision;
            
            obj.data.beta_est = obj.beta_cells(:, :, end);                  % collect final estimates
            obj.data.b_est = obj.beta_mean(end, :);
            obj.data.D_NTS = obj.D(:, :, 1);
            obj.data.D_GTS = obj.D(:, :, end);
                                                                            % predictions (cell-specific + population)
            obj.data.fitted2 = obj.system.integrate(obj.data.beta_est, obj.data);
            obj.data.population = obj.system.integrate(obj.data.b_est, obj.data);
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