classdef Optimization < handle   
    
    properties (Constant)
        options_quadprog = optimoptions('quadprog', 'Algorithm', 'active-set', 'Display', 'off', 'OptimalityTolerance', 1e-12);
        options_noisepar = optimoptions('fmincon', 'Display', 'off', 'OptimalityTolerance', 1e-3);
        options_initialization = optimoptions('fmincon', 'Display', 'off', 'StepTolerance', 1e-2);
    end
    

    methods (Static)
        function [opt, H] = QPGLS(design, response, covariance, x0, lb, ub, prior)
            if nargin < 7, prior = struct('mean', 0, 'prec', 0); end

            P = size(design, 2);                                            % number of parameters
            if nargin <= 5, lb = zeros(1, P); ub = Inf(1, P); end
            precision = svdinv(covariance);
            
            Y = reshape(response, [], 1);
            X = design;
            H = X' * precision * X + prior.prec;                            % quadratic programming problem
            H = (H+H')/2;                                                   % force symmetry
            f = -X' * precision * Y - prior.prec * prior.mean;
            if any(f ~= real(f), 'all') || any(f ~= real(f))
                disp(1)
            end                                                             % solve
            opt = quadprog(H, f, [], [], [], [], lb, ub, x0, Optimization.options_quadprog)';
        end
        
        
        function coef = noise_parameters(start, pred, measurements)
            ws = warning('off', 'MATLAB:nearlySingularMatrix');
            assert(numel(start) == 2)
            
            pred = flatten(pred);
            measurements = flatten(measurements);
            
            L = @(p) sum(log(p(1) + p(2).*pred.^2) + (measurements - pred).^2 ./ (p(1) + p(2).*pred.^2));
            coef = fmincon(L, start, [], [], [], [], [0 0], [Inf Inf], [], Optimization.options_noisepar);
            
            warning(ws)
        end


        function beta = initialize(SS, lb, ub, nstart)                  % Initial numerical optimization on population average
            value = Inf;
            for start = 1:nstart
                logb0 = rand(1, numel(lb)) .* (log(ub) - log(lb)) + log(lb);% uniformly sample in log parameter space
                [opt_new, value_new] = fmincon(@(logb0) SS(exp(logb0)), logb0, [], [], [], [], ...
                                               log(lb), log(ub), [], Optimization.options_initialization);
                if value_new < value, value = value_new; opt = opt_new; end
            end
            beta = exp(opt);                                                % transform back to normal scale
        end
    end
end


