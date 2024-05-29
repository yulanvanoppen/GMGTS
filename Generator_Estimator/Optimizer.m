classdef Optimizer < handle   
    
    properties (Constant)
        options_quadprog = optimoptions('quadprog', 'Algorithm', 'active-set', 'Display', 'off');
        options_newton = optimoptions('fmincon', 'Display', 'off', 'OptimalityTolerance', 1e-3);
%         options_newton = optimoptions('fmincon', 'Display', 'off', 'StepTolerance', 1e-2);
    end
    
    methods (Static)
        function beta = optimize(varargin)
            if nargin == 5
                beta = Optimizer.QPGLS(varargin{:});
            else
                beta = lsqnonneg(varargin{:})';
            end
        end
        

        function [opt, H] = QPGLS(design, response, covariance, x0, lb, ub, prior)
            if nargin < 7, prior = struct('mean', 0, 'prec', 0); end

            P = size(design, 2);                                    % number of parameters
            if nargin <= 5, lb = zeros(1, P); ub = Inf(1, P); end
            precision = svdinv(covariance);
            
            Y = reshape(response, [], 1);
            X = design;
            H = X' * precision * X + prior.prec;                    % quadratic programming problem
            H = (H+H')/2;                                           % force symmetry
            f = -X' * precision * Y - prior.prec * prior.mean;
            if any(f ~= real(f), 'all') || any(f ~= real(f))
                disp(1)
            end                                                     % solve
            opt = quadprog(H, f, [], [], [], [], lb, ub, x0, Optimizer.options_quadprog)';
        end
        
        
        function p = newton_raphson(p, x, y)
            ws = warning('off', 'MATLAB:nearlySingularMatrix');
            assert(numel(p) == 2)
            
            x = flatten(x);
            y = flatten(y);
            
            L = @(p) sum(log(p(1) + p(2).*x.^2) + (y - x).^2 ./ (p(1) + p(2).*x.^2));
            p = fmincon(L, p, [], [], [], [], [0 0], [Inf Inf], [], Optimizer.options_newton);
            
            warning(ws)
        end
    end
end


