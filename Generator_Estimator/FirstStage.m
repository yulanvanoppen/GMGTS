classdef FirstStage < handle
    properties (SetAccess = private)
        data                                                                % data aggregate struct
        system                                                              % nested object controlling the ODE system
        settings                                                            % hyperparameters and user input
        
        T                                                                   % number of time points
        L                                                                   % number of observed states
        N                                                                   % population size
        
        beta_fs                                                             % cell-specific parameter estimates
        V                                                                   % spline/ode gradient error variances
        varXdX                                                              % state and gradient covariances
        varbeta                                                             % uncertainty estimates of beta
        
        convergence_steps                                                   % relative iteration steps
        not_converged                                                       % indices of cells that have not yet converged
        
        smoothed_fitted                                                     % combination of predicted (hidden)
        dsmoothed_fitted                                                    % and smoothed (observed) states/gradients)
        fitted_fs                                                           % predicted states and gradients
        dfitted_fs
    end
    
    
    methods
        function obj = FirstStage(data, system, settings)               % Constructor
            obj.data = data;                                                % experimental and smoothing data
            obj.system = system;                                            % ODE system
            obj.settings = settings;                                        % hyperparameters
            
            [obj.T, obj.L, obj.N] = size(data.smoothed);                    % extract from trajectory dimensions
            
            obj.convergence_steps = ones(1, obj.N);                         % ensure no cells are considered converged
            obj.not_converged = 1:obj.N;
            
            obj.V = repmat(eye(system.K * obj.T), 1, 1, obj.N);
            obj.varXdX = zeros(2 * system.K * obj.T, 2 * system.K * obj.T, obj.N);
            obj.varbeta = repmat(eye(system.P), 1, 1, obj.N);
            
            if obj.L < system.K
                obj.initialize()                                            % numerically optimize
                obj.integrate(settings.niter == 0)                           % save solution if 'converged'
            end
        end
        
        
        function output = optimize(obj)                                 % Main optimization function
                                                                            % substitute smoothed measurements
            obj.smoothed_fitted(:, obj.data.observed, :) = obj.data.smoothed;
            obj.dsmoothed_fitted(:, obj.data.observed, :) = obj.data.dsmoothed;
            
            if obj.L == obj.system.K
                disp('Full observation GMGTS: Gradient matching only')
            else
                disp('Partial observation GMGTS: Extended iterative scheme')
            end

            if obj.settings.niter == 0, obj.estimate_covariances(1), end    % assume (premature) convergence
            
            for iter = 1:obj.settings.niter
                beta_old = obj.beta_fs;

                if iter > 1 && obj.L < obj.system.K, obj.integrate; end     % ODE integration
                obj.estimate_covariances(iter);                             % residual covariance estimation
                obj.update_parameters(iter);                                % gradient matching

                if obj.system.P > 1                                         % compute relative iteration steps
                    obj.convergence_steps = vecnorm((beta_old - obj.beta_fs)') ./ vecnorm(beta_old');
                else
                    obj.convergence_steps = abs(beta_old' - obj.beta_fs') ./ abs(beta_old');
                end                                                         % cells that have not converged 
                obj.not_converged = find(obj.convergence_steps >= obj.settings.tol);
                
                fprintf('%3d (%.1e, %.1e, %.1e)\n', iter, ...               % print iteration step quantiles .5, .9, 1.0
                        median(obj.convergence_steps), quantile(obj.convergence_steps, .9), max(obj.convergence_steps));
                
                                                                            % break when 90% of cells have converged
                if quantile(obj.convergence_steps, .9) < obj.settings.tol || iter == obj.settings.niter
                    if obj.L < obj.system.K                                 % estimate final uncertainty of beta
                        obj.estimate_covariances(iter, true);
                    else
                        obj.covariances_full_beta(iter);
                    end
                    break
                end
            end
            
            obj.beta_fs(obj.beta_fs < 1e-15) = 0;                           % polish final estimates
            obj.extract_estimates();                                        % save estimates and fitted trajectories
            output = obj.data;                                              % return data appended with FS results
        end
        
        
        function initialize(obj)                                        % Initial numerical optimization on population average
%             options = optimoptions('fmincon', 'Display', 'off', 'StepTolerance', 1e-2);
%             value = Inf;
%             for start = 1:obj.settings.nstart                               % uniformly sample in log parameter space
%                 logb0 = rand(1, obj.system.P) .* (log(obj.settings.ub) - log(obj.settings.lb)) + log(obj.settings.lb);
%                 [opt_new, value_new] = fmincon(@(logb0) obj.squares_sum(exp(logb0)), logb0, [], [], [], [], ...
%                                                log(obj.settings.lb), log(obj.settings.ub), [], options);
%                 if value_new < value, value = value_new; opt = opt_new; end
%             end
%             obj.beta_fs = repmat(exp(opt), obj.N, 1);                       % copy to each cell
            
            beta_init = Optimization.initialize(@obj.squares_sum, obj.settings.lb, obj.settings.ub, obj.settings.nstart);
            obj.beta_fs = repmat(beta_init, obj.N, 1);                       % copy to each cell
        end
            

        function ss = squares_sum(obj, b)                               % Weighted um of squared differences
            ss = Inf;
            try                                                             % compute fitted trajectories
                solution = obj.system.integrate(b, obj.data, obj.data.t_data);
                solution = solution(:, obj.data.observed, :);               % sum of squares on observed states
                ss = sum((solution - obj.data.traces).^2 ./ mean(obj.data.traces, [1 3]).^2, 'all');
                                                                            % add prior term (if any)
                ss = ss + (b' - obj.settings.prior.mean)' * (obj.settings.prior.prec * obj.settings.prior.mult) ...
                                                          * (b' - obj.settings.prior.mean);
            catch ME
                disp(ME)
            end
        end
            
        
        function update_parameters(obj, iter)                           % Update (cell-specific) parameter estimates
            for i = obj.not_converged                                       % model slopes from splines and integrations
                design = obj.system.g(obj.smoothed_fitted(2:end-1, :, i), obj.data.t(2:end-1));
                const = obj.system.const(obj.smoothed_fitted(2:end-1, :, i), obj.data.t(2:end-1));
                response = obj.dsmoothed_fitted(2:end-1, :, i) - const;     % spline and integration slopes
    
                                                                            % disregard interval ends
                nonzero_ind = reshape((2:obj.T-1)' + obj.T*(0:obj.system.K-1), 1, []);
                variances = obj.V(nonzero_ind, nonzero_ind, i);
                weights = reshape(sqrt(obj.settings.weights(2:end-1, :)), [], 1);
                variances = variances ./ weights ./ weights';
                
                if iter == 1 && obj.L == obj.system.K                       % omit starting point for full observation
                    initial = [];                                           % on the first iteration
                else
                    initial = obj.beta_fs(i, :);
                end                                                         % constrained GLS using quadratic programming
                obj.beta_fs(i, :) = Optimization.QPGLS(design, response, variances, initial, ...
                                                        obj.settings.lb, obj.settings.ub, obj.settings.prior);
            end
        end
        
                                        
        function estimate_covariances(obj, rep, converged)              % Update error covariance function parameters
            if nargin < 3, converged = false; end 
            
            if obj.L == obj.system.K
                obj.covariances_full;                                       % full observation
            else
                obj.covariances_partial(rep, converged);                    % partial observation
            end
        end
        
        
        function integrate(obj, converged)                              % Hidden states integration
            if nargin == 1, converged = false; end
            
            obj.fitted_fs = obj.system.integrate(obj.beta_fs, obj.data);    % integrate system and compute rhs
            obj.fitted_fs = max(1e-12, obj.fitted_fs);                      % force positive
            obj.dfitted_fs = obj.system.rhs(obj.fitted_fs, obj.data.t, obj.beta_fs);
            unobserved = ~ismember(1:obj.system.K, obj.data.observed);      % replace unobserved states by ODE integrations
            
                                                                            % substitute integrated hidden states and gradients
            obj.smoothed_fitted(:, unobserved, :) = obj.fitted_fs(:, unobserved, :); 
            obj.dsmoothed_fitted(:, unobserved, :) = obj.dfitted_fs(:, unobserved, :);
            
            if converged                                                    % store smooth versions upon convergence
                obj.data.fitted_fs_fine = max(1e-12, obj.system.integrate(obj.beta_fs, obj.data, obj.data.t_fine));
                obj.data.dfitted_fs_fine = obj.system.rhs(obj.data.fitted_fs_fine, obj.data.t_fine, obj.beta_fs);
            end
        end
        
        
        function extract_estimates(obj)                                 % Extract results 
            obj.integrate(true);                                            % compute fitted cell trajectories
            obj.data.beta_fs = obj.beta_fs;                                 % store results
            obj.data.fitted_fs = obj.fitted_fs;
            obj.data.dfitted_fs = obj.dfitted_fs;
            obj.data.smoothed_fitted_fs = obj.smoothed_fitted;
            obj.data.dsmoothed_fitted_fs = obj.dsmoothed_fitted;
            obj.data.V = obj.V;
            obj.data.varbeta = obj.varbeta;
            obj.data.convergence_steps = obj.convergence_steps;
            obj.data.converged = setdiff(1:obj.N, obj.not_converged);
        end
        
        
        function covariances_full(obj)                                  % Residual covariances for full observation
            Z = blkdiag(obj.data.basis{:})';                                % B-spline basis for all states
            Z_fs = blkdiag(obj.data.basis_fs{:})';                          % equivalent with first stage time point grid
            dZ_fs = blkdiag(obj.data.dbasis_fs{:})' / range(obj.data.t);    % corresponding differentiated basis
            
                                                                            % RHS Jacobian for each cell from smoothed measurements
            df_dX_all = obj.system.df(obj.data.smoothed, obj.data.t, obj.beta_fs);
            
            for i = obj.not_converged                                       % estimated measurement error variances
                S = diag(max(reshape(obj.data.variances_sm(:, :, i), 1, []), 1e-7));
                var_delta = svdinv(Z' * (S \ Z));                           % spline coefficient uncertainty
                var_smooth = [Z_fs; dZ_fs] * var_delta * [Z_fs' dZ_fs'];    % smoothed measurements covariance matrix
                
                df = zeros(obj.system.K * obj.T);                           % residual covariance approximation
                for j = 1:obj.T                                             % restructure RHS Jacobian
                    df(j + obj.T*(0:obj.system.K-1), j + obj.T*(0:obj.system.K-1)) = df_dX_all(j + obj.T*(0:obj.system.K-1), :, i);
                end
                
                residual_gradient = [-df eye(obj.L*obj.T)];                 % gradient for delta method
                
                regulator = max(1e-12, max(abs(obj.data.smoothed(:, :, i)) / 1e6, [], 1));
                regulator = reshape(repmat(regulator, obj.T, 1), 1, []);
                                                                            % regulated delta method approximation
                obj.V(:, :, i) = residual_gradient * var_smooth * residual_gradient' + diag(regulator);
            end
        end


        function covariances_full_beta(obj, rep)                        % Parameter uncertainties for full observation
            indices_t = 2:obj.T-1;                                          % time indices without interval ends
            indices_tk = reshape(indices_t' + obj.T*(0:obj.system.K-1), 1, []); % across states
            indices_2tk = reshape(indices_tk' + [0 obj.system.K*obj.T], 1, []); % across state and gradient components
                
            Z = blkdiag(obj.data.basis{:})';                                % B-spline basis for all states
            Z_fs = blkdiag(obj.data.basis_fs{:})';                          % equivalent with first stage time point grid
            dZ_fs = blkdiag(obj.data.dbasis_fs{:})' / range(obj.data.t);    % corresponding differentiated basis

            if rep > 1                                                      % compute g(.), h(.), and their Jacobians for each cell
                g_all = obj.system.g(obj.data.smoothed, obj.data.t);
                dg_all = obj.system.dg(obj.data.smoothed, obj.data.t);
                const_all = obj.system.const(obj.data.smoothed, obj.data.t);
                dconst_all = obj.system.dconst(obj.data.smoothed, obj.data.t);
            end
            
            for i = 1:obj.N                                                 % estimated measurement error variances
                S = diag(max(reshape(obj.data.variances_sm(:, :, i), 1, []), 1e-7));
                var_delta = svdinv(Z' * (S \ Z));                           % spline coefficient uncertainty
                var_smooth = [Z_fs; dZ_fs] * var_delta * [Z_fs' dZ_fs'];    % smoothed measurements covariance matrix
                
                dX = reshape(obj.data.dsmoothed(indices_t, :, i), [], 1);   % left-hand side
                H = reshape(const_all(indices_t, :, i), [], 1);             % constant part wrt parameters
                G = g_all(indices_tk, :, i);                                % linear part wrt parameters

                Vinv = zeros(obj.system.K * obj.T);                         % inverse residual covariance matrix estimate
                Vinv(indices_tk, indices_tk) = svdinv(obj.V(indices_tk, indices_tk, i));
                
                                                                            % construct d[beta]/d[dX] and d[beta]/d[X]
                                                                            % cf. section "S6 First-stage uncertainty estimates for fully observed systems"
                Th = G' * Vinv(indices_tk, indices_tk) * G + obj.settings.prior.prec;
                Thinv = svdinv(Th);
                dbeta_ddX = Thinv * G' * Vinv(indices_tk, indices_tk);      % partial with respect to dX
                
                Xi = G' * Vinv(indices_tk, indices_tk) * (dX - H) + obj.settings.prior.prec * obj.settings.prior.mean;
                Thinv_Xi = Thinv * Xi;
                Pi = zeros(obj.system.P, obj.T, obj.system.K);
                Psi = zeros(obj.system.P, obj.T, obj.system.K);
                for k = 1:obj.system.K
                    for j = 2:obj.T-1
                        dg_dxk = dg_all(:, :, j, k, i);
                        Vinv_j = Vinv(j + obj.T*(0:obj.system.K-1), indices_tk);

                        Pi(:, j, k) = (dg_dxk' * Vinv_j * G + G' * Vinv_j' * dg_dxk) * Thinv_Xi;
                        Psi(:, j, k) = dg_dxk' * Vinv_j * (dX-H) - G' * Vinv_j' * dconst_all(j + obj.T*(0:obj.system.K-1), k, i);
                    end
                end                                                         % partial with respect to X
                dbeta_dX = Thinv * reshape(Psi(:, indices_t, :) - Pi(:, indices_t, :), obj.system.P, []);

                                                                            % delta method approximation
                obj.varbeta(:, :, i) = [dbeta_dX dbeta_ddX] * var_smooth(indices_2tk, indices_2tk) * [dbeta_dX dbeta_ddX]';
            end
        end
        
        
        function covariances_partial(obj, rep, converged)               % Parameter uncertainties and residual covariances
            if nargin < 3, converged = false; end                       % for partial observation
            if converged, cells = 1:obj.N; else, cells = obj.not_converged; end

            hidden = setdiff(1:obj.system.K, obj.data.observed);            % hidden/observed state/time indices
            ind_hid = reshape((1:obj.T)'+(obj.T*(hidden-1)), [], 1);
            ind_obs = reshape((1:obj.T)'+(obj.T*(obj.data.observed-1)), [], 1);
            
            indices_t = 2:obj.T-1;                                          % time indices without interval ends
            indices_tk = reshape(indices_t' + obj.T*(0:obj.system.K-1), 1, []); % across states
            indices_2tk = reshape(indices_tk' + [0 obj.system.K*obj.T], 1, []); % across state and gradient components
            
            Z = blkdiag(obj.data.basis{:})';                                % B-spline basis for all states
            Z_fs = blkdiag(obj.data.basis_fs{:})';                          % equivalent with first stage time point grid
            dZ_fs = blkdiag(obj.data.dbasis_fs{:})' / range(obj.data.t);    % corresponding differentiated basis

            if rep > 1                                                      % compute g(.), h(.), and their Jacobians for each cell
                g_all = obj.system.g(obj.smoothed_fitted, obj.data.t);
                dg_all = zeros(obj.system.K, obj.system.P, obj.T, obj.system.K, obj.N);
                dg_all(:, :, :, :, cells) = obj.system.dg(obj.smoothed_fitted(:, :, cells), obj.data.t);
                const_all = obj.system.const(obj.smoothed_fitted, obj.data.t);
                dconst_all = obj.system.dconst(obj.smoothed_fitted, obj.data.t);
            end
                                                                            % RHS Jacobian for each cell from smoothed measurements
            df_dX_all = obj.system.df(obj.smoothed_fitted, obj.data.t, obj.beta_fs);

            beta_permuted = permute(obj.beta_fs, [3 2 1]);                  % vectorized finite difference approximations
            epsilon = max(1e-8, beta_permuted * .001);                      % of solution and RHS parameter sensitivities
            beta_pm_eps = beta_permuted + epsilon .* kron([1; -1], eye(obj.system.P));

            traces_pm_eps = zeros(obj.T, obj.system.K, obj.system.P, 2, obj.N);
            for i = cells
                for p = 1:obj.system.P                                      % numerically integrate at parameter component pm eps
                    traces_pm_eps(:, :, p, 1, i) = obj.system.integrate(beta_pm_eps(p, :, i), obj.data, obj.data.t, 1e-4);
                    traces_pm_eps(:, :, p, 2, i) = obj.system.integrate(beta_pm_eps(p+end/2, :, i), obj.data, obj.data.t, 1e-4);
                end
            end
                                                                            % restructure for rhs evaluations
            beta_pm_eps = reshape(permute(beta_pm_eps, [2 1 3]), obj.system.P, [])';
            gradients_pm_eps = reshape(obj.system.rhs(reshape(traces_pm_eps, obj.T, obj.system.K, []), obj.data.t, ...
                                                    beta_pm_eps), obj.T, obj.system.K, obj.system.P, 2, []);

            eps_denom = 1./reshape(2*epsilon, 1, 1, obj.system.P, 1, obj.N);% finite difference approximations
            dF_dbeta_all = permute((traces_pm_eps(:, :, :, 1, :) - traces_pm_eps(:, :, :, 2, :)) .* eps_denom, [1:3 5 4]);
            df_dbeta_all = permute((gradients_pm_eps(:, :, :, 1, :) - gradients_pm_eps(:, :, :, 2, :)) .* eps_denom, [1:3 5 4]);
            
            for i = cells
                dF_dbeta = reshape(dF_dbeta_all(:, :, :, i), [], obj.system.P);
                df_dbeta = reshape(df_dbeta_all(:, :, :, i), [], obj.system.P);

                if rep == 1                                                 % smoothing covariance
                    S = diag(max(reshape(obj.data.variances_sm(:, :, i), 1, []), 1e-7));
                    var_delta = svdinv(Z' * (S \ Z));
                    var_XOdXO = [Z_fs; dZ_fs] * var_delta * [Z_fs; dZ_fs]';

                    obj.varbeta(:, :, i) = .5^2 * diag(obj.beta_fs(i, :).^2);   % initial (large) parameter uncertainty
                    var_XdXint = [dF_dbeta; df_dbeta] * obj.varbeta(:, :, i) * [dF_dbeta; df_dbeta]';

                    var_XdX = zeros(2 * obj.system.K * obj.T);              % full state/gradient covariance
                    var_XdX(ind_obs + [0 end/2], ind_obs + [0 end/2]) = var_XOdXO;
                    var_XdX(ind_hid + [0 end/2], ind_hid + [0 end/2]) = var_XdXint(ind_hid + [0 end/2], ...
                                                                                       ind_hid + [0 end/2]);
                    obj.varXdX(indices_2tk, indices_2tk, i) = var_XdX(indices_2tk, indices_2tk);

                else
                                                                            % left-hand side
                    dX = reshape(obj.dsmoothed_fitted(indices_t, :, i), [], 1);
                    H = reshape(const_all(indices_t, :, i), [], 1);         % constant part wrt parameters
                    G = g_all(indices_tk, :, i);                            % linear part wrt parameters
                    
                    Vinv = zeros(obj.system.K * obj.T);                     % inverse residual covariance matrix estimate
                    Vinv(indices_tk, indices_tk) = svdinv(obj.V(indices_tk, indices_tk, i));
                    
                                                                            % construct d[beta]/d[dX] and d[beta]/d[X]
                                                                            % cf. section "S6 First-stage uncertainty estimates for fully observed systems"
                    Th = G' * Vinv(indices_tk, indices_tk) * G + obj.settings.prior.prec;
                    Thinv = svdinv(Th);
                    dbeta_ddX = Thinv * G' * Vinv(indices_tk, indices_tk);  % partial with respect to dX

                    Xi = G' * Vinv(indices_tk, indices_tk) * (dX - H) + obj.settings.prior.prec * obj.settings.prior.mean;
                    Thinv_Xi = Thinv * Xi;
                    Pi = zeros(obj.system.P, obj.T, obj.system.K);
                    Psi = zeros(obj.system.P, obj.T, obj.system.K);
                    for k = 1:obj.system.K
                        for j = 2:obj.T-1
                            dg_dxk = dg_all(:, :, j, k, i);
                            Vinv_j = Vinv(j + obj.T*(0:obj.system.K-1), indices_tk);

                            Pi(:, j, k) = (dg_dxk' * Vinv_j * G + G' * Vinv_j' * dg_dxk) * Thinv_Xi;
                            Psi(:, j, k) = dg_dxk' * Vinv_j * (dX-H) - G' * Vinv_j' * dconst_all(j + obj.T*(0:obj.system.K-1), k, i);
                        end
                    end                                                     % partial with respect to X
                    dbeta_dX = Thinv * reshape(Psi(:, indices_t, :) - Pi(:, indices_t, :), obj.system.P, []);
                    
                                                                            % delta method approximation
                    obj.varbeta(:, :, i) = [dbeta_dX dbeta_ddX] * obj.varXdX(indices_2tk, indices_2tk, i) * [dbeta_dX dbeta_ddX]';
                    if converged, return, end

                                                                            % iterative update of state/gradient covariance
                    identity_O_indices = zeros(obj.L * obj.T, obj.system.K * obj.T);
                    identity_O_indices(:, ind_obs) = eye(obj.L * obj.T);
                    delta_XdX = zeros(2 * obj.system.K * obj.T);
                    
                    delta_XdX(ind_hid,         indices_2tk) = dF_dbeta(ind_hid, :)  * [dbeta_dX dbeta_ddX];
                    delta_XdX(ind_hid + end/2, indices_2tk) = df_dbeta(ind_hid, :) * [dbeta_dX dbeta_ddX];

                    delta_XdX(ind_obs,         :) = [identity_O_indices   zeros(obj.L * obj.T, obj.system.K * obj.T)];
                    delta_XdX(ind_obs + end/2, :) = [zeros(obj.L * obj.T, obj.system.K * obj.T)   identity_O_indices];
                    
                    obj.varXdX(indices_2tk, indices_2tk, i) = delta_XdX(indices_2tk, indices_2tk) * obj.varXdX(indices_2tk, indices_2tk, i) ...
                                                                                                  * delta_XdX(indices_2tk, indices_2tk)';
                end
                
                df_dX = zeros(obj.system.K * obj.T);                        % residual covariance approximation
                for j = 1:obj.T                                             % restructure RHS Jacobian
                    df_dX(j + obj.T*(0:obj.system.K-1), j + obj.T*(0:obj.system.K-1)) = df_dX_all(j + obj.T*(0:obj.system.K-1), :, i);
                end
                residual_gradient = [-df_dX eye(obj.system.K*obj.T)];       % delta method approximation
                V_unregulated = residual_gradient(indices_tk, indices_2tk) * obj.varXdX(indices_2tk, indices_2tk, i) ...
                                                                           * residual_gradient(indices_tk, indices_2tk)';
                
                regulator = max(1e-12, max(abs(obj.smoothed_fitted(:, :, i)) / 1e6, [], 1));
                regulator = reshape(repmat(regulator, obj.T-2, 1), 1, []);
                                                                            % regulated delta method approximation
                obj.V(indices_tk, indices_tk, i) = V_unregulated + diag(regulator);
            end
        end
    end
end








