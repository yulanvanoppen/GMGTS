classdef Smoother < handle
    properties (Access = private)
        T                                                                   % number of time points
        T_fs                                                                % number of optimized time points
        L                                                                   % system dimension
        N                                                                   % number of cells
    end
    
    properties (SetAccess = private)
        data                                                                % general data and output container
        system
        settings
        bsplines                                                             % nested object controlling splines from B-splines
        
        iteration
        delta
        B                                                                   % B-spline basis evaluated at t
        dB                                                                  % corresponding first derivative
        ddB                                                                 % corresponding second derivative
        B_fine
        dB_fine
        smoothed                                                            % smoothed data
        dsmoothed                                                           % smoothed derivative
        variances_sm
        variances_fs
        sigma
        tau
    end
    
    methods
        function obj = Smoother(data, system, settings)                 % Constructor
            obj.data = data;                                                % store experiment data
            obj.system = system;
            obj.settings = settings;
            [obj.T, obj.L, obj.N] = size(obj.data.traces);                  % extract trajectory dimensions
            
            if settings.autoknots, obj.place_knots(); end
            
            obj.bsplines = cell(1, obj.L);
            obj.delta = cell(1, obj.L);
            for k = 1:obj.L                                             % instantiate spline
                obj.bsplines{k} = BSpline(obj.settings.order, obj.settings.knots{k});
                obj.delta{k} = zeros(obj.bsplines{k}.card, obj.N);
            end
            
            obj.variances_sm = zeros(obj.T, obj.L, obj.N);
            obj.variances_fs = zeros(length(settings.grid), obj.L, obj.N);  % measurement error variances
            obj.sigma = zeros(1, obj.L);
            obj.tau = zeros(1, obj.L);
        end

        function update_knots(obj, knots)
            if iscell(knots), obj.settings.knots = knots; end
            for k = 1:obj.L
                if ~iscell(knots), obj.settings.knots{k} = knots; end
                obj.bsplines{k} = BSpline(obj.settings.order, obj.settings.knots{k});
                obj.delta{k} = zeros(obj.bsplines{k}.card, obj.N);
            end
        end
        
        
        function output = smooth(obj)                                   % Smooth trajectory data
            if obj.settings.interactive
                app = smoothing_new(obj);
                waitfor(app.FinishButton, 'UserData')
                if isvalid(app)
                    obj.data = app.smoother.data;
                    delete(app)
                else
                    obj.optimize()
                end
            else
                obj.optimize()                                              % return extended data struct
            end
            
            obj.data.t_data = obj.data.t;
            obj.data.T_data = obj.data.T;
            obj.data.t = obj.settings.grid;
            obj.data.T = length(obj.data.t);
            
            output = obj.data;
        end
        
        
        function optimize(obj)
            tic
            
            grid = (obj.data.t - obj.data.t(1)) / range(obj.data.t);        % scaled time grid
            grid_fine = linspace(0, 1, 201);                                % scaled optimized time grid
            grid_fs = (obj.settings.grid - obj.data.t(1)) / range(obj.data.t);
            obj.T_fs = length(grid_fs);
            
            obj.B = cell(1, obj.L);
            obj.dB = cell(1, obj.L);
            obj.B_fine = cell(1, obj.L);
            obj.dB_fine = cell(1, obj.L);
            obj.data.basis_fs = cell(1, obj.L);
            obj.data.dbasis_fs = cell(1, obj.L);
            
            for k = 1:obj.L
                obj.B{k} = obj.bsplines{k}.basis(grid);                     % B-spline basis evaluated at t
                obj.dB{k} = obj.bsplines{k}.basis_diff(grid);               % corresponding first derivative
                obj.data.basis_fs{k} = obj.bsplines{k}.basis(grid_fs);      % save for later use
                obj.data.dbasis_fs{k} = obj.bsplines{k}.basis_diff(grid_fs); 
                obj.B_fine{k} = obj.bsplines{k}.basis(grid_fine);           % save for smooth plotting
                obj.dB_fine{k} = obj.bsplines{k}.basis_diff(grid_fine);
            end

            obj.data.basis = obj.B;                                         % save for later use
            obj.data.dbasis = obj.dB;

            for iter = 1:obj.settings.niter                          % FWLS spline coefficient estimates
                obj.iteration = iter;
                delta_old = obj.delta; 
                obj.update_coefficients();                     
                obj.update_variances();
                if eucl_rel(delta_old', obj.delta') < obj.settings.tol      % check convergence
                    break
                end
            end
            
            obj.data.variances_sm = obj.variances_sm;
            obj.data.variances = obj.variances_fs;                          % save variance estimates
            obj.fit();                                                      % compute splines
            
            obj.data.toc_sm = toc;
        end
        
                                            
        function update_coefficients(obj)                               % Spline coefficients from current variances
            for i = 1:obj.N                                                 % W: correcting weights from variances
                for k = 1:obj.L                                             % penalty: lambda * L2(spline basis/coefficients)
                    W = diag(max(obj.variances_sm(:, k, i), 1e-7).^-1);     % delta_ik: LS estimate
                    obj.delta{k}(:, i) = svdinv(obj.B{k} * W * obj.B{k}') * obj.B{k} * W * obj.data.traces(:, k, i);
                end
            end
        end
        
        
        function update_variances(obj)                                  % Estimate variances from current splines
            for k = 1:obj.L
                predicted = obj.B{k}' * reshape(obj.delta{k}, obj.bsplines{k}.card, obj.N);
                predicted_fs = obj.data.basis_fs{k}' * reshape(obj.delta{k}, obj.bsplines{k}.card, obj.N);
                design = [ones(obj.N * obj.T, 1) flatten(predicted).^2];
                design_fs = [ones(obj.N * obj.T_fs, 1) flatten(predicted_fs).^2];
                response = (flatten(predicted) - flatten(obj.data.traces(:, k, :))).^2;
                
                coefficients = Optimizer.optimize(design, response);        % find nonzero LS estimates for theta
                if sum(coefficients) == 0, coefficients(1) = mean(response); end
                                                                            % optimize further iteratively
                coefficients = Optimizer.newton_raphson(coefficients, predicted, obj.data.traces(:, k, :));
                if sum(coefficients) == 0, coefficients(1) = mean(response); end
                
                obj.sigma(k) = coefficients(1);                                % store opt and compute variances
                obj.tau(k) = coefficients(2);
                obj.variances_sm(:, k, :) = reshape(design * coefficients', obj.T, 1, obj.N);
                obj.variances_fs(:, k, :) = reshape(design_fs * coefficients', obj.T_fs, 1, obj.N);
            end
        end
        
        
        function fit(obj)                                               % Fitted splines from estimated coefficients
            obj.data.smoothed = zeros(obj.T_fs, obj.L, obj.N);              % smoothed data
            obj.data.dsmoothed = zeros(obj.T_fs, obj.L, obj.N);             % smoothed derivative
            obj.data.smoothed_fine = zeros(201, obj.L, obj.N);              % same for finer time grid
            obj.data.dsmoothed_fine = zeros(201, obj.L, obj.N);             % smoothed derivative
            
            for i = 1:obj.N                                                 % coefficients to spline interpolants
                for k = 1:obj.L
                    obj.data.smoothed(:, k, i) = obj.data.basis_fs{k}' * obj.delta{k}(:, i);
                    obj.data.dsmoothed(:, k, i) = obj.data.dbasis_fs{k}' * obj.delta{k}(:, i) / range(obj.data.t);
                    obj.data.smoothed_fine(:, k, i) = obj.B_fine{k}' * obj.delta{k}(:, i);
                    obj.data.dsmoothed_fine(:, k, i) = obj.dB_fine{k}' * obj.delta{k}(:, i) / range(obj.data.t);
                end
            end
%             a = linspace(1, 0, size(obj.data.smoothed, 1))';
%             obj.data.smoothed = a .* obj.data.smoothed + (1-a) .* obj.data.original(:, 2, :);
%             obj.data.dsmoothed = a .* obj.data.dsmoothed + (1-a) .* obj.data.doriginal(:, 2, :);
            
            obj.data.smoothed = max(obj.data.smoothed, 1e-12);
        end
        
        
        function place_knots(obj)
            t = obj.data.t';

            base = false(obj.T-2, obj.L);                   % 3 EQUIDISTANT BASE KNOTS
            [~, base_ind] = min(abs(t - linspace(t(1), t(end), 5)));
            base_ind = unique(base_ind);
            base(base_ind(2:end-1), :) = true; 

            d21 = t(2:end-1) - t(1:end-2);                  % FINITE DIFFERENCE DERIVATIVE APPROXIMATIONS
            d31 = t(3:end) - t(1:end-2);                        % with unequal time step
            d32 = t(3:end) - t(2:end-1);
            y1 = obj.data.traces(1:end-2, :, :);
            y2 = obj.data.traces(2:end-1, :, :);
            y3 = obj.data.traces(3:end, :, :);

            dy = (y3 - y2) ./ (d32 + d21);                      % unequally spaced finite differences
            ddy = 2 * (y1./d21./d31 - y2./d32./d21 + y3./d32./d31);

            zscore_dy = mean(dy, 3) ./ std(dy, 0, 3);           % means and standard deviations across cells
            zscore_ddy = mean(ddy, 3) ./ std(ddy, 0, 3);
                                                            % ADD ESTIMATED PEAKS AND TROUGHS
            crossings = logical(abs(diff(sign(zscore_dy))) == 2);
            crossings(end+1, :) = false;
            for state = 1:obj.L                                 % find approximate points where dy = 0
                for idx = 1:size(crossings, 1)-1
                    if crossings(idx, state)                    % select dy before or after crossing closest to zero
                        subset_abs_mdy = abs(zscore_dy(idx:idx+1, state));
                        crossings(idx:idx+1, state) = subset_abs_mdy == min(subset_abs_mdy);
                    end
                end
            end
            peaks_troughs = crossings & (abs(zscore_ddy) > .25);% filter noise crossings

            base_peaks_troughs = base;                          % add peaks/troughs moving base points if needed
            for state = 1:obj.L
                for idx = 1:size(peaks_troughs, 1)
                    if peaks_troughs(idx, state)                % remove adjacent base points
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
                                                            % ADD END OF INITIAL FAST DYNAMICS
            curvature = abs(zscore_ddy) > 2*std(zscore_ddy);    % find first point where curvature settles
            fast_dynamics_end = ones(1, obj.L);
            for state = 1:obj.L                                 % first point with two or more regular curvatures
                n_fast_dynamics = find(diff([find(curvature(:, state)); obj.T+1]) > 2, 1);
                fast_dynamics_ind = find(curvature(:, state), n_fast_dynamics);
                fast_dynamics_end(state) = fast_dynamics_ind(end);
            end

            base_curvature = base_peaks_troughs;                % add fast dynamics end moving base points if needed
            for state = 1:obj.L
                idx = fast_dynamics_end(state);                 % remove adjacent base points
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

            subset = true(obj.T, obj.L);                    % ADD KNOTS AT INTERVAL ENDS
            subset(2:end-1, :) = base_curvature;
            for state = 1:obj.L                                 % normalize knots to [0, 1]
                obj.settings.knots{state} = (t(subset(:, state))' - t(1)) / range(t);
            end
        end
    end
end