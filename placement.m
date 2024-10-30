function knots = placement(data)
    t = data.t';
    T = data.T;
    L = length(data.observed);
    
    % BASE KNOTS
    base = false(data.T-2, L);
    [~, base_ind] = unique(min(abs(data.t' - linspace(data.t(1), data.t(end), 5))));
    base(base_ind(2:end-1), :) = true; 
    
    % FINITE DIFFERENCE DERIVATIVE APPROXIMATIONS
    d21 = t(2:end-1) - t(1:end-2);                      
    d31 = t(3:end) - t(1:end-2);                        % with unequal time step
    d32 = t(3:end) - t(2:end-1);
    y1 = data.traces(1:end-2, :, :);
    y2 = data.traces(2:end-1, :, :);
    y3 = data.traces(3:end, :, :);
    
    dy = (y3 - y2) ./ (d32 + d21);                      % unequally spaced finite differences
    ddy = 2 * (y1./d21./d31 - y2./d32./d21 + y3./d32./d31);
    
                                                        % means and standard deviations across cells
    my_norm = mean(data.traces, 3) ./ std(data.traces, 0, [1 3]);
    mdy_norm = mean(dy, 3) ./ std(dy, 0, 3);
    mddy_norm = mean(ddy, 3) ./ std(ddy, 0, 3);

    curvature_threshold = 2*std(mddy_norm);
    
    figure(1)
    tiledlayout(1, 3)
    nexttile, plot(t, my_norm), legend, grid
    nexttile, plot(t(2:end-1), mdy_norm), legend, grid
    nexttile, plot(t(2:end-1), mddy_norm), legend, grid
    
    % CRITERIA: PEAKS, TROUGHS, AND HIGH CURVATURE
    crossings = logical(abs(diff(sign(mdy_norm))) == 2);
    crossings(end+1, :) = false;
    for state = 1:L
        for idx = 1:size(crossings, 1)-1
            if crossings(idx, state)
                subset_abs_mdy = abs(mdy_norm(idx:idx+1, state));
                crossings(idx:idx+1, state) = subset_abs_mdy == min(subset_abs_mdy);
            end
        end
    end
    peaks_troughs = crossings & (abs(mddy_norm) > .25);
    
    base_peaks_troughs = base, 1, L;
    for state = 1:L
        for idx = 1:size(peaks_troughs, 1)
            if peaks_troughs(idx, state)
                if idx > 1 && idx < data.T-2 && base(idx-1, state) && base(idx+1, state)
                    base_peaks_troughs(idx-1:idx+1, state) = [false true false];
                elseif idx > 1 && base(idx-1, state)
                    base_peaks_troughs(idx-1:idx, state) = [false true];
                elseif idx < data.T-2 && base(idx+1, state)
                    base_peaks_troughs(idx:idx+1, state) = [true false];
                else
                    base_peaks_troughs(idx, state) = true;
                end
            end
        end
    end
    
    curvature = abs(mddy_norm) > curvature_threshold;
    fast_dynamics_end = ones(1, L);
    for state = 1:L
        n_fast_dynamics = find(diff([find(curvature(:, state)); data.T+1]) > 2, 1);
        fast_dynamics_ind = find(curvature(:, state), n_fast_dynamics);
        fast_dynamics_end(state) = fast_dynamics_ind(end);
    end
    
    base_curvature = base_peaks_troughs;
    for state = 1:L
        idx = fast_dynamics_end(state);
        if idx > 1 && idx < data.T-2 && base_peaks_troughs(idx-1, state) && base_peaks_troughs(idx+1, state)
            base_curvature(idx-1:idx+1, state) = [false true false];
        elseif idx > 1 && base_peaks_troughs(idx-1, state)
            base_curvature(idx-1:idx, state) = [false true];
        elseif idx < data.T-2 && base_peaks_troughs(idx+1, state)
            base_curvature(idx:idx+1, state) = [true false];
        else
            base_curvature(idx, state) = true;
        end
    end
    
    crossings2 = logical(abs(diff(sign(mddy_norm))) == 2);
    crossings2(end+1, :) = false;
    for state = 1:L
        for idx = 1:size(crossings2, 1)-1
            if crossings2(idx, state)
                subset_abs_mddy = abs(mddy_norm(idx:idx+1, state));
                crossings(idx:idx+1, state) = subset_abs_mddy == min(subset_abs_mddy);
            end
        end
    end
    
%     inflections = false(size(crossings2));
%     for state = 1:L
%         prev = 0;
%         idx = find(crossings2(:, state), 1);
%         while idx < size(crossings2, 1)
%             next = idx + find(crossings2(idx+1:end, state), 1);
%             if isempty(next), next = size(crossings2, 1); end
%             prev_idx = mddy_norm(prev+1:idx, state);
%             idx_next = mddy_norm(idx+1:next, state);
%             [~, max_prev] = max(abs(prev_idx));
%             [~, max_next] = max(abs(idx_next));
%             before_crossing = mddy_norm(prev+max_prev, state);
%             after_crossing = mddy_norm(idx+max_next, state);
%             if ~isempty(before_crossing) && ~isempty(after_crossing) ...
%                && (max(abs(before_crossing), abs(after_crossing)) > std(mddy_norm(:, state))/4) ...
%                && sign(before_crossing) ~= sign(after_crossing)
%                inflections(idx, state) = true;
%             end
%             prev = idx;
%             idx = next;
%         end
%     end
%     for state = 1:L
%         idx = 1;
%         while idx < size(inflections, 1)
%             if inflections(idx, state) && ~diff(inflections(idx:idx+1, state))
%                 subset_abs_mddy = abs(mddy_norm(idx:idx+1, state));
%                 inflections(idx:idx+1, state) = subset_abs_mddy == min(subset_abs_mddy);
%                 idx = 1;
%                 continue
%             end
%             idx = idx+1;
%         end
%     end
%     inflections
    
%     % FILTER DOUBLE KNOTS
%     combined = peaks_troughs | curvature | inflections;
%     for state = 1:L
%         idx = 1;
%         while idx < size(mddy_norm, 1)
%             if combined(idx, state) && ~diff(combined(idx:idx+1, state))
%                 subset_abs_mddy = abs(mddy_norm(idx:idx+1, state));
%                 combined(idx:idx+1, state) = subset_abs_mddy == max(subset_abs_mddy);
%                 idx = 1;
%                 continue
%             end
%             idx = idx+1;
%         end
%     end

%     for state = 1:L
%         state
%         idx = 1;
%         while idx < size(mddy_norm, 1)
%             if combined(idx, state)
%                 next_false = idx + find(~combined(idx+1:end, state), 1);
%                 if isempty(next_false), next_false = size(mddy_norm, 1) + 1; end
%                 if next_false - idx > 1
%                     subset_abs_mddy = abs(mddy_norm(idx:next_false-1, state));
%                     combined(idx:next_false-1, state) = subset_abs_mddy == max(subset_abs_mddy);
%                     idx = next_false;
%                     continue
%                 end
%             end
%             idx = idx+1;
%         end
%     end
    
%     % PENALIZED INTERVAL BASED ON HIGH INITIAL CURVATURE
%     penalized = zeros(2, L);
%     for state = 1:L
%         first_high_curvature = find(curvature(:, state), 1);
%         if first_high_curvature < data.T/4
%             slow_dynamics = first_high_curvature + find(~curvature(first_high_curvature+1:end, state), 1);
%             penalized(:, state) = [t(slow_dynamics+1) t(end)];
%         else
%             penalized(:, state) = [t(1) t(end)];
%         end
%     end
    
    % REMOVED KNOTS AT SECOND TO LAST TIME POINT
    subset = true(T, L);
    subset(2:end-1, :) = base_curvature;
    knots = cell(1, L);
    for state = 1:L
        knots{state} = t(subset(:, state))';
    end
end