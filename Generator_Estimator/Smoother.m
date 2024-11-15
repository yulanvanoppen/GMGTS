classdef Smoother < handle
    properties (SetAccess = private)
        data                                                                % general data and output container
        system                                                              
        settings                                                            % hyperparameters and user input
        bsplines                                                            % nested objects controlling spline smoothing
        
        T                                                                   % number of data time points
        T_fs                                                                % number of first stage time points
        T_fine                                                              % number of time points for smooth plots
        L                                                                   % system dimension
        N                                                                   % number of cells
        
        scaled                                                              % data time grid scaled to [0, 1]
        scaled_fs                                                           % first stage scaled time grid 
        scaled_fine                                                         % fine time scaled time grid for plotting
        
        B                                                                   % B-spline basis evaluated at data time grid
        dB                                                                  % corresponding first derivative
        B_fine                                                              % equivalent basis and its derivative 
        dB_fine                                                             % at fine time grid for plotting

        delta                                                               % B-spline coefficients

        sigma                                                               % additive measurement error variance 
        tau                                                                 % multiplicative measurement error variance 

        variances_sm                                                        % measurement error variances at data time grid
        variances_fs                                                        % and first stage time grid
    end
    
    
    methods
        function obj = Smoother(data, system, settings)                 % Constructor
            obj.data = data;                                                % store data, ODE system, and hyperparameters
            obj.system = system;                                            
            obj.settings = settings;
            [obj.T, obj.L, obj.N] = size(obj.data.traces);                  % extract trajectory dimensions
            obj.T_fs = length(settings.t_fs);

            obj.scaled = (obj.data.t - obj.data.t(1)) / range(obj.data.t);  % scale to [0, 1]
            obj.scaled_fs = (obj.settings.t_fs - obj.data.t(1)) / range(obj.data.t);
            obj.scaled_fine = (obj.data.t_fine - obj.data.t(1)) / range(obj.data.t);
            
            if settings.autoknots, obj.place_knots(); end
            
            obj.bsplines = cell(1, obj.L);                                  % initialize splines and coefficients
            obj.delta = cell(1, obj.L);
            obj.update_knots(obj.settings.knots);
            
            obj.sigma = zeros(1, obj.L);                                    % initialize measurement error variances
            obj.tau = zeros(1, obj.L);
            obj.variances_sm = zeros(obj.T, obj.L, obj.N);
            obj.variances_fs = zeros(length(settings.t_fs), obj.L, obj.N);
        end

        
        function update_knots(obj, knots, states)
            if nargin < 3, states = 1:obj.L; end
            for k = states                                                  % copy knots, reset bases and coefficients
                obj.settings.knots{k} = knots{k};                           
                obj.bsplines{k} = BSpline(obj.settings.order, obj.settings.knots{k});
                obj.delta{k} = zeros(obj.bsplines{k}.card, obj.N);
            end
            
            [obj.B, obj.dB, obj.B_fine, obj.dB_fine, ...
                    obj.data.basis_fs, obj.data.dbasis_fs] = deal(cell(1, obj.L));
            
            for k = 1:obj.L
                obj.B{k} = obj.bsplines{k}.basis(obj.scaled);               % B-spline basis evaluated at t
                obj.dB{k} = obj.bsplines{k}.basis_diff(obj.scaled);         % corresponding first derivative
                obj.data.basis_fs{k} = obj.bsplines{k}.basis(obj.scaled_fs);% save for later use
                obj.data.dbasis_fs{k} = obj.bsplines{k}.basis_diff(obj.scaled_fs); 
                obj.B_fine{k} = obj.bsplines{k}.basis(obj.scaled_fine);     % save for smooth plotting
                obj.dB_fine{k} = obj.bsplines{k}.basis_diff(obj.scaled_fine);
            end

            obj.data.basis = obj.B;                                         % save for later use
            obj.data.dbasis = obj.dB;
        end
        
        
        function output = smooth(obj)                                   % Smooth trajectory data
            if obj.settings.interactive                                     % if interactive smoothing, open app
                app = smoothing_new(obj);
                waitfor(app.FinishButton, 'UserData')                       % wait until finished
                if isvalid(app), delete(app), else, obj.optimize(), end     % close app
            else
                obj.optimize()                                              % estimate coefficients and noise parameters
            end
            
            obj.data.t_data = obj.data.t;                                   % prime time point data for first stage
            obj.data.T_data = obj.data.T;
            obj.data.t = obj.settings.t_fs;
            obj.data.T = length(obj.data.t);
            
            output = obj.data;                                              % return data extended with smoothing results
        end
        
        
        function optimize(obj, states)
            if nargin == 1, states = 1:obj.L; end                           % don't overwrite computing time when 
            if length(states) == obj.L, tic, end                            % the app smoothes single states
            
            for iter = 1:obj.settings.niter                                 % FWLS spline coefficient estimates
                delta_old = obj.delta; 
                obj.update_coefficients(states);                            % WLS estimate
                obj.update_variances(states);                               % noise parameter MLE
                if eucl_rel(delta_old', obj.delta') < obj.settings.tol      % check convergence
                    break
                end
            end
            
            obj.data.variances_sm = obj.variances_sm;                       % save variance estimates wrt
            obj.data.variances = obj.variances_fs;                          % measurement and first-stage time points
            obj.fit();                                                      % compute splines
            
            if length(states) == obj.L, obj.data.toc_sm = toc; end
        end
        
                                            
        function update_coefficients(obj, states)                       % Spline coefficients from current variances
            for i = 1:obj.N                                                 % W: correcting weights from variances
                for k = states                                              % penalty: lambda * L2(spline basis/coefficients)
                    W = diag(max(obj.variances_sm(:, k, i), 1e-7).^-1);     % delta_ik: LS estimate
                    obj.delta{k}(:, i) = svdinv(obj.B{k} * W * obj.B{k}') * obj.B{k} * W * obj.data.traces(:, k, i);
                end
            end
        end
        
        
        function update_variances(obj, states)                          % Estimate variances from current splines
            for k = states                                                  % smoothed measurements
                predicted = obj.B{k}' * reshape(obj.delta{k}, obj.bsplines{k}.card, obj.N);
                design = [ones(obj.N*obj.T, 1) flatten(predicted).^2];      % columns for additive and multiplicative noise
                response = (flatten(predicted) - flatten(obj.data.traces(:, k, :))).^2;
                
                coefficients = lsqnonneg(design, response)';                 % initialize nonzero LS estimates for noise parameters
                if sum(coefficients) == 0, coefficients(1) = mean(response); end
                                                                            % optimize further iteratively
                coefficients = Optimization.noise_parameters(coefficients, predicted, obj.data.traces(:, k, :));
                if sum(coefficients) == 0, coefficients(1) = mean(response); end
                
                obj.sigma(k) = coefficients(1);                             % store optimum and compute variances
                obj.tau(k) = coefficients(2);
                obj.variances_sm(:, k, :) = reshape(design * coefficients', obj.T, 1, obj.N);
                
                                                                            % equivalent for first stage time grid
                predicted_fs = obj.data.basis_fs{k}' * reshape(obj.delta{k}, obj.bsplines{k}.card, obj.N);
                design_fs = [ones(obj.N * obj.T_fs, 1) flatten(predicted_fs).^2];
                obj.variances_fs(:, k, :) = reshape(design_fs * coefficients', obj.T_fs, 1, obj.N);
            end
        end
        
        
        function fit(obj)                                               % Fitted splines from estimated coefficients
            obj.data.smoothed = zeros(obj.T_fs, obj.L, obj.N);              % smoothed measurements
            obj.data.dsmoothed = zeros(obj.T_fs, obj.L, obj.N);             % corresponding gradients
            obj.data.smoothed_fine = zeros(81, obj.L, obj.N);               % equivalent on fine time grid for plotting
            obj.data.dsmoothed_fine = zeros(81, obj.L, obj.N);
            
            for i = 1:obj.N
                for k = 1:obj.L                                             % multiply coefficients with bases
                    obj.data.smoothed(:, k, i) = obj.data.basis_fs{k}' * obj.delta{k}(:, i);
                    obj.data.dsmoothed(:, k, i) = obj.data.dbasis_fs{k}' * obj.delta{k}(:, i) / range(obj.data.t);
                    obj.data.smoothed_fine(:, k, i) = obj.B_fine{k}' * obj.delta{k}(:, i);
                    obj.data.dsmoothed_fine(:, k, i) = obj.dB_fine{k}' * obj.delta{k}(:, i) / range(obj.data.t);
                end
            end
            
            obj.data.smoothed = max(obj.data.smoothed, 1e-12);              % clip nonpositive values
        end
        
        
        function place_knots(obj)                                       % Knot placement heuristic
            t = obj.data.t';

            base = false(obj.T-2, obj.L);                                   % start with three equidistant interior knots
            [~, base_ind] = min(abs(t - linspace(t(1), t(end), 5)));        % choose closest data time points
            base_ind = unique(base_ind);                                    % remove any duplicates
            base(base_ind(2:end-1), :) = true;                              % disregard interval ends

            d21 = t(2:end-1) - t(1:end-2);                                  % finite difference derivative approximations
            d31 = t(3:end) - t(1:end-2);                                    % with unequal time step
            d32 = t(3:end) - t(2:end-1);
            y1 = obj.data.traces(1:end-2, :, :);
            y2 = obj.data.traces(2:end-1, :, :);
            y3 = obj.data.traces(3:end, :, :);
            dy = (y3 - y2) ./ (d32 + d21);
            ddy = 2 * (y1./d21./d31 - y2./d32./d21 + y3./d32./d31);

            zscore_dy = mean(dy, 3) ./ std(dy, 0, 3);                       % normalized mean slopes and curvature estimates
            zscore_ddy = mean(ddy, 3) ./ std(ddy, 0, 3);
            
            crossings = logical(abs(diff(sign(zscore_dy))) == 2);           % add estimated peaks and troughs
            crossings(end+1, :) = false;
            for state = 1:obj.L                                             % find approximate points where dy = 0
                for idx = 1:size(crossings, 1)-1
                    if crossings(idx, state)                                % select dy closest to zero near crossing
                        subset_abs_mdy = abs(zscore_dy(idx:idx+1, state));
                        crossings(idx:idx+1, state) = subset_abs_mdy == min(subset_abs_mdy);
                    end
                end
            end
            peaks_troughs = crossings & (abs(zscore_ddy) > .25);            % filter crossings due to noise

            base_peaks_troughs = base;                                      % add peaks/troughs moving base points if needed
            for state = 1:obj.L
                for idx = 1:size(peaks_troughs, 1)
                    if peaks_troughs(idx, state)                            % remove adjacent base points
                        if idx > 1 && idx < obj.T-2 && base(idx-1, state) && base(idx+1, state)
                            base_peaks_troughs(idx-1:idx+1, state) = [false true false];
                        elseif idx > 1 && base(idx-1, state)
                            base_peaks_troughs(idx-1:idx, state) = [false true];
                        elseif idx < obj.T-2 && base(idx+1, state)
                            base_peaks_troughs(idx:idx+1, state) = [true false];
                        else
                            base_peaks_troughs(idx, state) = true;
                        end
                    end
                end
            end
                                                                            
            curvature = abs(zscore_ddy) > 2*std(zscore_ddy);                % add end of initial fast dynamics (if any)
            fast_dynamics_end = zeros(1, obj.L);
            for state = 1:obj.L                                             % first point with two or more regular curvatures
                if any(curvature(:, state))
                    n_fast_dynamics = find(diff([find(curvature(:, state)); obj.T+2]) > 2, 1);
                    fast_dynamics_ind = find(curvature(:, state), n_fast_dynamics);
                    fast_dynamics_end(state) = fast_dynamics_ind(end);
                end
            end

            base_curvature = base_peaks_troughs;                            % add fast dynamics end moving existing points if needed
            for state = 1:obj.L
                idx = fast_dynamics_end(state);                             % remove adjacent existing points
                if idx == 0, continue, end
                if idx > 1 && idx < obj.T-2 && base_peaks_troughs(idx-1, state) && base_peaks_troughs(idx+1, state)
                    base_curvature(idx-1:idx+1, state) = [false true false];
                elseif idx > 1 && base_peaks_troughs(idx-1, state)
                    base_curvature(idx-1:idx, state) = [false true];
                elseif idx < obj.T-2 && base_peaks_troughs(idx+1, state)
                    base_curvature(idx:idx+1, state) = [true false];
                else
                    base_curvature(idx, state) = true;
                end
            end

            subset = true(obj.T, obj.L);                                    % add knots at interval ends
            subset(2:end-1, :) = base_curvature;
            for state = 1:obj.L                                             % normalize knots to [0, 1]
                obj.settings.knots{state} = (t(subset(:, state))' - t(1)) / range(t);
            end
        end
    end
end