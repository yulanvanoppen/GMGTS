classdef SecondStage < handle
    properties (SetAccess = immutable)
        init                                    % initial values
        t                                       % time grid
        dt_int                                  % integration time step
        dt_dat                                  % data time step
        varying                                 % indices of parameters that may vary among cells
        
        nrep                                    % maximum number of iterations
        tol                                     % desired error tolerance
        
        T                                       % number of time points
        L                                       % #observed states
        K                                       % system dimension
        N                                       % population size
        P                                       % parameter vector length
    end
    
    properties (SetAccess = private)
        system                                  % nested object controlling the ODE system
        settings
        
        converged_fs
        precision                               % asymptotic precision estimate
        beta_cells                              % cell-specific parameter iterations
        beta_mean                               % population parameter iterations
        D                                       % random effect covariance iterations
    end
    
    properties (SetAccess = public)
        data
        
        b_est                                   % final cell-specific parameter estimates
        beta_est                                % final population parameter estimates
        
        D_NTS                                   % naive random effect covariance estimate
        D_GTS                                   % final random effect covariance estimate
        
        fitted2                                 % fitted cell trajectories
        population                              % fitted population trajectories
    end
    
    methods                                 % Constructor
        function obj = SecondStage(data, system, settings)
            obj.data = data;
            obj.system = system;
            obj.settings = settings;
            
            obj.init = data.init;
            obj.t = data.t;
            obj.varying = data.varying;
            
            obj.nrep = settings.nrep;
            obj.tol = settings.tol;
            
            obj.converged_fs = data.convergence_steps < .5;
            
                                                % extract from trajectory dimensions
            [obj.T, obj.L, obj.N] = size(obj.data.traces);
            obj.K = system.K;
            obj.P = system.P;
            
            obj.beta_cells = data.beta_fs;   % generate initial estimate from first stage results
%             obj.beta_cells = data.beta_fs(obj.converged_fs, :);
            obj.beta_mean = mean(data.beta_fs);
            obj.D = cov(data.beta_fs(:, obj.varying));
        end
        
                                            % Main optimization function
        function output = optimize(obj)         % max. absolute distance among argument elements
            obj.estimate_precisions_new();
            
            for rep = 1:obj.nrep                % current D^-1
%                 fprintf('%3d ', rep)
%                 if mod(rep, 10) == 0, fprintf('\n'), end
                
                certain = diag(obj.D(:, :, rep)) < 1e-7;
                Dinv_rep = zeros(size(obj.D(:, :, rep)));
                Dinv_rep(~certain, ~certain) = svdinv(obj.D(~certain, ~certain, rep));
                Dinv_rep(certain, certain) = diag(Inf(1, sum(certain)));
                
                                                % iterate until convergence:
                obj.E_step(rep, Dinv_rep)       %   - refine cell-specific estimates
                obj.M_step(rep, Dinv_rep)       %   - refine population estimates and covariances
                
                                                % break upon convergence
                if eucl_rel(obj.beta_mean(rep, :), obj.beta_mean(rep+1, :)) < obj.tol && ...
                   eucl_rel(obj.D(:, :, rep), obj.D(:, :, rep+1), true) < obj.tol
                    break
                end
            end
            
%             fprintf('\n')
            
            obj.extract_estimates();            % store final estimates publically
            output = obj.data;
        end
        
                                            % Expectation step of EM algorithm
        function E_step(obj, rep, Dinv_rep)     % abbreviate for clarity
            beta_rep = obj.beta_mean(rep, obj.varying)';

                                                % copy to retain non-varying parameters
            obj.beta_cells(:, :, rep+1) = obj.beta_cells(:, :, rep);

            for i = 1:obj.N                     % cell i's initial beta and Ci^-1
                init_i = obj.beta_cells(i, obj.varying, 1)';
                Cinv_i = obj.precision(:, :, i);
                
                certain = isinf(diag(Cinv_i)) | isinf(diag(Dinv_rep));
                
                                                % update estimate
                obj.beta_cells(i, obj.varying(certain), rep+1) = init_i(certain);
                obj.beta_cells(i, obj.varying(~certain), rep+1) = ((Cinv_i(~certain, ~certain) ...
                                                                    + Dinv_rep(~certain, ~certain)) ...
                                                                 \ (Cinv_i(~certain, ~certain) * init_i(~certain) ...
                                                                    + Dinv_rep(~certain, ~certain) * beta_rep(~certain)))';
            end
        end

                                            % Maximization step of EM algorithm
        function M_step(obj, rep, Dinv_rep)
            warning('off','MATLAB:singularMatrix')
                                                % update population parameters and abbreviate
            obj.beta_mean(rep+1, :) = obj.beta_mean(rep, :);
            obj.beta_mean(rep+1, :) = mean(obj.beta_cells(:, :, rep+1));
            beta_rep = obj.beta_mean(rep+1, obj.varying)';

                                                % collect terms of D
            summands = zeros(length(obj.varying), length(obj.varying), obj.N);

            for i = 1:obj.N                     % abbreviate cell i's (current) beta and C^-1
                beta_i = obj.beta_cells(i, obj.varying, rep+1)';
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


            obj.data.D_NTS = obj.D(:, :, 1);                % - random effect covariance matrices
            obj.data.D_GTS = obj.D(:, :, end);
                                               % compute final fitted cell trajectories
            obj.data.fitted2 = obj.system.integrate(obj.data.beta_est, obj.data, obj.data.t_fine);
            
                                               % compute population mean trajectory
            obj.data.population = obj.system.integrate(obj.data.b_est, obj.data, obj.data.t_fine);
        end
        
        
        function estimate_precisions_new(obj)     % inverse covariance matrix estimates
%             if obj.L == obj.system.K
%                 sqweights = sqrt(obj.settings.weights); 
%                 W = diag(reshape(sqweights, [], 1));                      
%                 g_all = obj.system.g(obj.data.smoothed_fitted_fs, obj.data.t);
%                 obj.precision = zeros(length(obj.data.varying), length(obj.data.varying), obj.N);
%             end
            
            for i = 1:obj.N
%                 if obj.L < obj.system.K
                    covariance = obj.data.variances_beta_fs(:, :, i);
%                 else
%                     V = obj.data.variances_fs(:, :, i);
%                     Vinv = svdinv(V);
%                     G = g_all(:, obj.data.varying, i);
%                     J = G' * W * Vinv * W * G;
%                     covariance = svdinv(J);
%                 end
                
                certain = diag(covariance) < max(diag(covariance)) / 1e6;
                obj.precision(~certain, ~certain, i) = svdinv(covariance(~certain, ~certain));
                obj.precision(certain, certain, i) = diag(Inf(1, sum(certain)));
                
%                 if i == 1 && any(diag(covariance) > 1)
%                     disp(1)
%                 end
            end
        end
    end
end