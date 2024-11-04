classdef FirstStage < handle
    properties (SetAccess = private)
        data                                                                % data aggregate struct
        system                                                              % nested object controlling the ODE system
        settings                                                            % weights and loop controls
        
        T                                                                   % number of time points
        L                                                                   % number of observed states
        N                                                                   % population size
        
        beta_fs                                                             % cell-specific parameter estimates
        V                                                                   % spline/ode gradient error variances
        varXdX                                                              % state and gradient covariances
        varbeta                                                             % uncertainty estimates of beta
        
        convergence_steps                                                   % relative iteration steps
        not_converged                                                       % indices of cells that have not yet converged
        
        smoothed_fitted                                                     % combination of predictions (hidden states)
        dsmoothed_fitted                                                    % and smoothing (observed states) + derivatives
        fitted_fs                                                           % predictions + derivatives
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
                obj.integrate(settings.nrep == 0)                           % save solution if 'converged'
            end
        end
        
        
        function output = optimize(obj)                                 % Main optimization function
            obj.smoothed_fitted(:, obj.data.observed, :) = obj.data.smoothed;
            obj.dsmoothed_fitted(:, obj.data.observed, :) = obj.data.dsmoothed;
            
            if obj.L == obj.system.K                                        % full observation
                disp('Full observation GMGTS: Gradient matching only')
            else
                disp('Partial observation GMGTS: Extended iterative scheme')
            end

            if obj.settings.nrep == 0, obj.estimate_covariances(1), end
            
            for rep = 1:obj.settings.nrep
                beta_old = obj.beta_fs;

                if rep > 1 && obj.L < obj.system.K, obj.integrate; end      % ODE integration
                obj.estimate_covariances(rep);                              % residual covariance estimation
                obj.update_parameters(rep);                                 % gradient matching

                if obj.system.P > 1                                         % compute relative iteration steps
                    obj.convergence_steps = vecnorm((beta_old - obj.beta_fs)') ./ vecnorm(beta_old');
                else
                    obj.convergence_steps = abs(beta_old' - obj.beta_fs') ./ abs(beta_old');
                end
                obj.not_converged = find(obj.convergence_steps >= obj.settings.tol);
                
                fprintf('%3d (%.1e, %.1e, %.1e)\n', rep, ...
                        median(obj.convergence_steps), quantile(obj.convergence_steps, .9), max(obj.convergence_steps));
                
                if quantile(obj.convergence_steps, .9) < obj.settings.tol || rep == obj.settings.nrep
                    if obj.L < obj.system.K
                        obj.estimate_covariances(rep, true);
                    else
                        obj.covariances_full_beta(rep);
                    end
                    obj.integrate(true);
                    break
                end
            end
            
            obj.beta_fs(obj.beta_fs < 1e-15) = 0;                           % polish final estimates
            obj.extract_estimates();                                        % save estimates and fitted trajectories
            output = obj.data;                                              % return expanded data container
        end
        
        
        function initialize(obj)
            options = optimoptions('fmincon', 'Display', 'off', 'StepTolerance', 1e-2);
%             options = optimoptions('fmincon', 'Display', 'off', 'OptimalityTolerance', 1e-2);
            value = Inf;
            for start = 1:obj.settings.nstart                               % multistart on first run
                logb0 = rand(1, obj.system.P) .* (log(obj.settings.ub) - log(obj.settings.lb)) + log(obj.settings.lb);
                [opt_new, value_new] = fmincon(@(logb0) obj.squares_sum(exp(logb0)), logb0, [], [], [], [], ...
                                               log(obj.settings.lb), log(obj.settings.ub), [], options);
                if value_new < value, value = value_new; opt = opt_new; end
            end
            disp("Initial beta:")
            disp(exp(opt))
            obj.beta_fs = repmat(exp(opt), obj.N, 1);
        end
            

        function ss = squares_sum(obj, b)              % Weighted um of squared differences
            ss = Inf;
            try
                solution = obj.system.integrate(b, obj.data, obj.data.t_data);
                solution = solution(:, obj.data.observed, :);
                ss = sum((solution - obj.data.traces).^2 ./ mean(obj.data.traces, [1 3]).^2, 'all');

                ss = ss + (b' - obj.settings.prior.mean)' * (obj.settings.prior.prec * obj.settings.prior.mult) ...
                                                          * (b' - obj.settings.prior.mean);
            catch ME
                disp(ME)
            end
        end
            
        
        function update_parameters(obj, rep)                        % Update (cell-specific) parameter estimates
            for i = obj.not_converged                                       % model slopes from splines + integrations
                design = obj.system.g(obj.smoothed_fitted(2:end-1, :, i), obj.data.t(2:end-1));
                const = obj.system.const(obj.smoothed_fitted(2:end-1, :, i), obj.data.t(2:end-1));
                response = obj.dsmoothed_fitted(2:end-1, :, i) - const;     % spline + integration slopes

                nonzero_ind = reshape((2:obj.T-1)' + obj.T*(0:obj.system.K-1), 1, []);
                variances = obj.V(nonzero_ind, nonzero_ind, i);
                weights = reshape(sqrt(obj.settings.weights(2:end-1, :)), [], 1);
                variances = variances ./ weights ./ weights';
                
                if rep == 1 && obj.L == obj.system.K
                    initial = [];
                else
                    initial = obj.beta_fs(i, :);
                end 
                obj.beta_fs(i, :) = Optimizer.QPGLS(design, response, variances, initial, ...
                                                    obj.settings.lb, obj.settings.ub, obj.settings.prior);
            end
        end
        
                                        
        function estimate_covariances(obj, rep, converged)              % Update error covariance function parameters
            if nargin < 3, converged = false; end 
            
            if obj.L == obj.system.K
                obj.covariances_full;
            else
                obj.covariances_partial(rep, converged);
            end
        end
        
        
        function integrate(obj, converged)                              % Hidden states integration
            if nargin == 1, converged = false; end
            
            obj.fitted_fs = obj.system.integrate(obj.beta_fs, obj.data);    % integrate system
            obj.fitted_fs = max(1e-12, obj.fitted_fs);
            obj.dfitted_fs = obj.system.rhs(obj.fitted_fs, obj.data.t, obj.beta_fs);    % compute corresponding rhs
            unobserved = ~ismember(1:obj.system.K, obj.data.observed);      % replace unobserved states by ODE integrations
            
            obj.smoothed_fitted(:, unobserved, :) = obj.fitted_fs(:, unobserved, :); 
            obj.dsmoothed_fitted(:, unobserved, :) = obj.dfitted_fs(:, unobserved, :);
            
            if converged
                obj.data.t_fine = linspace(obj.data.t(1), obj.data.t(end), 201);
                obj.data.fitted_fs_fine = max(1e-12, obj.system.integrate(obj.beta_fs, obj.data, obj.data.t_fine));
                obj.data.dfitted_fs_fine = obj.system.rhs(obj.data.fitted_fs_fine, obj.data.t_fine, obj.beta_fs);
            end
        end
        
        
        function extract_estimates(obj)                                 % Extract results 
            obj.integrate();                                                % compute fitted cell trajectories
            obj.data.beta_fs = obj.beta_fs;                                 % store results
            obj.data.fitted_fs = obj.fitted_fs;
            obj.data.dfitted_fs = obj.dfitted_fs;
            obj.data.smoothed_fitted_fs = obj.smoothed_fitted;
            obj.data.dsmoothed_fitted_fs = obj.dsmoothed_fitted;
            obj.data.V = obj.V;
            obj.data.varbeta = obj.varbeta;
            obj.data.convergence_steps = obj.convergence_steps;
        end
        
        
        function covariances_full(obj)
            Z = blkdiag(obj.data.basis{:})';
            Z_fs = blkdiag(obj.data.basis_fs{:})';
            dZ_fs = blkdiag(obj.data.dbasis_fs{:})' / range(obj.data.t);
            
            dRHS_all = obj.system.df(obj.data.smoothed, obj.data.t, obj.beta_fs);
            
            for i = obj.not_converged
                S = diag(reshape(obj.data.variances_sm(:, :, i) + 1e-12, 1, []));
                var_delta = svdinv(Z' * (S \ Z));
                var_smooth = [Z_fs; dZ_fs] * var_delta * [Z_fs' dZ_fs'];
                
                dF = zeros(obj.system.K * obj.T);
                for j = 1:obj.T
                    dF(j + obj.T*(0:obj.system.K-1), j + obj.T*(0:obj.system.K-1)) = dRHS_all(j + obj.T*(0:obj.system.K-1), :, i);
                end
                
                delta_h = [-dF eye(obj.L*obj.T)];
                
                regulator = max(1e-12, max(abs(obj.data.smoothed(:, :, i)) / 1e6, [], 1));
                regulator = reshape(repmat(regulator, obj.T, 1), 1, []);
                
                obj.V(:, :, i) = delta_h * var_smooth * delta_h' + diag(regulator);
            end
        end


        function covariances_full_beta(obj, rep)
            wt_ind = 2:obj.T-1;
            wtk_ind = reshape(wt_ind' + obj.T*(0:obj.system.K-1), 1, []);
            wtk_ind2 = reshape(wtk_ind' + [0 obj.system.K*obj.T], 1, []);
            
            states = obj.smoothed_fitted;
            dstates = obj.dsmoothed_fitted;
                
            Z = blkdiag(obj.data.basis{:})';
            Z_fs = blkdiag(obj.data.basis_fs{:})';
            dZ_fs = blkdiag(obj.data.dbasis_fs{:})' / range(obj.data.t);

            if rep > 1
                g_all = obj.system.g(states, obj.data.t);
                dg_all = obj.system.dg(states, obj.data.t);
                const_all = obj.system.const(states, obj.data.t);
                dconst_all = obj.system.dconst(states, obj.data.t);
            end
            
            for i = 1:obj.N
                S = diag(reshape(obj.data.variances_sm(:, :, i) + 1e-12, 1, []));
                var_delta = svdinv(Z' * (S \ Z));
                var_XdX = [Z_fs; dZ_fs] * var_delta * [Z_fs' dZ_fs'];
                
                dX = reshape(dstates(wt_ind, :, i), [], 1);
                const = reshape(const_all(wt_ind, :, i), [], 1);
                G = g_all(wtk_ind, :, i);

                Vinv = zeros(obj.system.K * obj.T);
                Vinv(wtk_ind, wtk_ind) = svdinv(obj.V(wtk_ind, wtk_ind, i));
                
                Th = G' * Vinv(wtk_ind, wtk_ind) * G + obj.settings.prior.prec;
                Thinv = svdinv(Th);
                dbeta_ddX = Thinv * G' * Vinv(wtk_ind, wtk_ind);
                
                Xi = G' * Vinv(wtk_ind, wtk_ind) * (dX - const) + obj.settings.prior.prec * obj.settings.prior.mean;
                Thinv_Xi = Thinv * Xi;
                Pi = zeros(obj.system.P, obj.T, obj.system.K);
                Psi = zeros(obj.system.P, obj.T, obj.system.K);
                for k = 1:obj.system.K
                    for j = 2:obj.T-1
                        dg_dxk = dg_all(:, :, j, k, i);
                        Vinv_j = Vinv(j + obj.T*(0:obj.system.K-1), wtk_ind);

                        Pi(:, j, k) = (dg_dxk' * Vinv_j * G + G' * Vinv_j' * dg_dxk) * Thinv_Xi;
                        Psi(:, j, k) = dg_dxk' * Vinv_j * (dX-const) - G' * Vinv_j' * dconst_all(j + obj.T*(0:obj.system.K-1), k, i);
                    end
                end
                dbeta_dX = Thinv * reshape(Psi(:, wt_ind, :) - Pi(:, wt_ind, :), obj.system.P, []);

                obj.varbeta(:, :, i) = [dbeta_dX dbeta_ddX] * var_XdX(wtk_ind2, wtk_ind2) * [dbeta_dX dbeta_ddX]';
            end
        end
        
        
        function covariances_partial(obj, rep, converged)
            if nargin < 3, converged = false; end
            if converged, cells = 1:obj.N; else, cells = obj.not_converged; end
            
            observed = obj.data.observed;
            hidden = 1:obj.system.K;
            hidden = hidden(~ismember(1:obj.system.K, observed));
            H_indices = reshape((1:obj.T)'+(obj.T*(hidden-1)), [], 1);
            O_indices = reshape((1:obj.T)'+(obj.T*(observed-1)), [], 1);
            
            wt_ind = 2:obj.T-1;
            wtk_ind = reshape(wt_ind' + obj.T*(0:obj.system.K-1), 1, []);
            wtk_ind2 = reshape(wtk_ind' + [0 obj.system.K*obj.T], 1, []);
%             zero_wt_ind = reshape(([1 obj.T])' + obj.T*(0:obj.system.K-1), 1, []);

            states = obj.smoothed_fitted;
            dstates = obj.dsmoothed_fitted;
            
            Z = blkdiag(obj.data.basis{:})';
            Z_fs = blkdiag(obj.data.basis_fs{:})';
            dZ_fs = blkdiag(obj.data.dbasis_fs{:})' / range(obj.data.t); 

            if rep > 1
                g_all = obj.system.g(states, obj.data.t);
                dg_all = zeros(obj.system.K, obj.system.P, obj.T, obj.system.K, obj.N);
                dg_all(:, :, :, :, cells) = obj.system.dg(states(:, :, cells), obj.data.t);
                const_all = obj.system.const(states, obj.data.t);
                dconst_all = obj.system.dconst(states, obj.data.t);
            end
            
            dRHS_all = obj.system.df(states, obj.data.t, obj.beta_fs);

            beta_transf = permute(obj.beta_fs, [3 2 1]);
            epsilon = max(1e-8, beta_transf * .001);
            beta_pm_eps = beta_transf + epsilon .* kron([1; -1], eye(obj.system.P));
            traces_pm_eps = zeros(obj.T, obj.system.K, obj.system.P, 2, obj.N);
            
            for i = cells
                for p = 1:obj.system.P
                    traces_pm_eps(:, :, p, 1, i) = obj.system.integrate(beta_pm_eps(p, :, i), obj.data, obj.data.t, 1e-4);
                    traces_pm_eps(:, :, p, 2, i) = obj.system.integrate(beta_pm_eps(p+end/2, :, i), obj.data, obj.data.t, 1e-4);
                end
            end
            
            beta_pm_eps = reshape(permute(beta_pm_eps, [2 1 3]), obj.system.P, [])';
            dtraces_pm_eps = reshape(obj.system.rhs(reshape(traces_pm_eps, obj.T, obj.system.K, []), obj.data.t, ...
                                                    beta_pm_eps), obj.T, obj.system.K, obj.system.P, 2, []);
            eps_denom = reshape(2*epsilon, 1, 1, obj.system.P, 1, obj.N).^-1;
            gradient_all = permute((traces_pm_eps(:, :, :, 1, :) - traces_pm_eps(:, :, :, 2, :)) .* eps_denom, [1:3 5 4]);
            dgradient_all = permute((dtraces_pm_eps(:, :, :, 1, :) - dtraces_pm_eps(:, :, :, 2, :)) .* eps_denom, [1:3 5 4]);
            
            for i = cells
                % DEFINITION OF dF_dbeta & ddF_dbeta
                dF_dbeta = reshape(gradient_all(:, :, :, i), [], obj.system.P);
                ddF_dbeta = reshape(dgradient_all(:, :, :, i), [], obj.system.P);

                if rep == 1
                    % DEFINITION OF Var XO
%                     S = diag(reshape(obj.data.variances_sm(:, :, i) + 1e-12, 1, []));
                    S = diag(max(reshape(obj.data.variances_sm(:, :, i), 1, []), 1e-7));
                    var_delta = svdinv(Z' * (S \ Z));
                    var_XOdXO = [Z_fs; dZ_fs] * var_delta * [Z_fs; dZ_fs]';

                    % DEFINITION OF Var beta
                    var_beta = .5^2 * diag(obj.beta_fs_init(1, :).^2);

                    obj.varbeta(:, :, i) = var_beta;
                    var_XdXint = [dF_dbeta; ddF_dbeta] * var_beta * [dF_dbeta; ddF_dbeta]';

                    % DEFINTION OF Var X
                    var_XdX = zeros(2 * obj.system.K * obj.T);
                    var_XdX(O_indices + [0 end/2], O_indices + [0 end/2]) = var_XOdXO;
                    var_XdX(H_indices + [0 end/2], H_indices + [0 end/2]) = var_XdXint(H_indices + [0 end/2], ...
                                                                                       H_indices + [0 end/2]);
                    obj.varXdX(wtk_ind2, wtk_ind2, i) = var_XdX(wtk_ind2, wtk_ind2);
                else
                    % DEFINITION OF dbeta_dXO and dbeta_ddXO
                    dX = reshape(dstates(wt_ind, :, i), [], 1);
                    const = reshape(const_all(wt_ind, :, i), [], 1);
                    G = g_all(wtk_ind, :, i);
                    
                    Vinv = zeros(obj.system.K * obj.T);
                    Vinv(wtk_ind, wtk_ind) = svdinv(obj.V(wtk_ind, wtk_ind, i));

%                     weights = reshape(sqrt(obj.settings.weights), [], 1);
%                     weights = reshape(sqrt(obj.settings.weights .* obj.settings.state_weights), [], 1);
%                     Vinv = svdinv(obj.V(:, :, i)) .* weights .* weights';
                    
%                     Th = G' * Vinv * G;
                    Th = G' * Vinv(wtk_ind, wtk_ind) * G + obj.settings.prior.prec;
                    Thinv = svdinv(Th);
                    dbeta_ddX = Thinv * G' * Vinv(wtk_ind, wtk_ind);
%                     dbeta_ddX(:, zero_wt_ind) = 0;

%                     Xi = G' * Vinv * (dX - const);
                    Xi = G' * Vinv(wtk_ind, wtk_ind) * (dX - const) + obj.settings.prior.prec * obj.settings.prior.mean;
                    Thinv_Xi = Thinv * Xi;
                    Pi = zeros(obj.system.P, obj.T, obj.system.K);
                    Psi = zeros(obj.system.P, obj.T, obj.system.K);
                    for k = 1:obj.system.K
                        for j = 2:obj.T-1
                            dg_dxk = dg_all(:, :, j, k, i);
                            Vinv_j = Vinv(j + obj.T*(0:obj.system.K-1), wtk_ind);

                            Pi(:, j, k) = (dg_dxk' * Vinv_j * G + G' * Vinv_j' * dg_dxk) * Thinv_Xi;
                            Psi(:, j, k) = dg_dxk' * Vinv_j * (dX-const) - G' * Vinv_j' * dconst_all(j + obj.T*(0:obj.system.K-1), k, i);
                        end
                    end
                    dbeta_dX = Thinv * reshape(Psi(:, wt_ind, :) - Pi(:, wt_ind, :), obj.system.P, []);
%                     dbeta_dX(:, zero_wt_ind) = 0;
                    
                    obj.varbeta(:, :, i) = [dbeta_dX dbeta_ddX] * obj.varXdX(wtk_ind2, wtk_ind2, i) * [dbeta_dX dbeta_ddX]';
                    if converged, return, end

                    % DEFINITION OF var_X
                    identity_O_indices = zeros(obj.L * obj.T, obj.system.K * obj.T);
                    identity_O_indices(:, O_indices) = eye(obj.L * obj.T);
                    delta_XdX = zeros(2 * obj.system.K * obj.T);
                    
                    delta_XdX(H_indices,         wtk_ind2) = dF_dbeta(H_indices, :)  * [dbeta_dX dbeta_ddX];
                    delta_XdX(H_indices + end/2, wtk_ind2) = ddF_dbeta(H_indices, :) * [dbeta_dX dbeta_ddX];

                    delta_XdX(O_indices,         :) = [identity_O_indices   zeros(obj.L * obj.T, obj.system.K * obj.T)];
                    delta_XdX(O_indices + end/2, :) = [zeros(obj.L * obj.T, obj.system.K * obj.T)   identity_O_indices];
                    
                    obj.varXdX(wtk_ind2, wtk_ind2, i) = delta_XdX(wtk_ind2, wtk_ind2) * obj.varXdX(wtk_ind2, wtk_ind2, i) * delta_XdX(wtk_ind2, wtk_ind2)';
                end
                
                % VARIANCE OF X, dX, AND dX - G(X)beta - H(X)
                dRHS = zeros(obj.system.K * obj.T);
                for j = 1:obj.T
                    dRHS(j + obj.T*(0:obj.system.K-1), j + obj.T*(0:obj.system.K-1)) = dRHS_all(j + obj.T*(0:obj.system.K-1), :, i);
                end
                delta_h = [-dRHS eye(obj.system.K * obj.T)];
                
                V_unregulated = delta_h(wtk_ind, wtk_ind2) * obj.varXdX(wtk_ind2, wtk_ind2, i) * delta_h(wtk_ind, wtk_ind2)';
                
                regulator = max(1e-12, max(abs(obj.smoothed_fitted(:, :, i)) / 1e6, [], 1));
                regulator = reshape(repmat(regulator, obj.T-2, 1), 1, []);
                obj.V(wtk_ind, wtk_ind, i) = V_unregulated + diag(regulator);
            end
        end
    end
end








