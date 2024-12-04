function [b, D] = EM_lognormal(beta, varbeta)
    ws = warning('error', 'MATLAB:nearlySingularMatrix');

    beta_lnorm = zeros(size(beta));
    varbeta_lnorm = zeros(size(varbeta));

    for i = 1:size(beta, 1)                                                 % lognormal approximation
        beta_i = beta(i, :)';
        varbeta_i = varbeta(:, :, i);
                                                                            % moment matching
%         obj.data.varbeta_lnorm(:, :, i) = log(1 + varbeta_i ./ (beta_i*beta_i'));
%         obj.data.beta_lnorm(i, :) = log(beta_i') - .5 * diag(obj.data.varbeta_lnorm(:, :, i))';

                                                                            % unbiased log-normal approximation
        diag_i = diag(log(.5 + sqrt(.25 + diag(diag(varbeta_i)./beta_i.^2))));
        full_i = log(1 + varbeta_i ./ (beta_i * beta_i') .* exp(-.5 * (diag_i + diag_i')));
        varbeta_lnorm(:, :, i) = full_i;
        beta_lnorm(i, :) = log(beta_i');
    end
    
    precision_lnorm = invert_uncertainties(varbeta_lnorm);
    
    [b, D] = EM(beta_lnorm, precision_lnorm);
    
    warning(ws);
 end


function [b, D] = EM(beta_lnorm, precision_lnorm)
    N = size(beta_lnorm, 1);

    b = mean(beta_lnorm);
    D = cov(beta_lnorm);
    beta_EM = beta_lnorm;
    
    for iter = 1:20
        b_old = b;
        D_old = D;
        Dinv = inv(D);
                                                                            % iterate until convergence:
        for i = 1:N                                                         % E-step
            beta_i = beta_lnorm(i, :)';
            Cinv_i = precision_lnorm(:, :, i);
            beta_EM(i, :) = ((Cinv_i+Dinv) \ (Cinv_i*beta_i + Dinv*b'))';   % update cell-specific parameters
        end
                                                                            % M-step
        b = mean(beta_EM);                                                  % update b
        summands = zeros(size(precision_lnorm));                            % collect terms of D
        for i = 1:N
            beta_i = beta_EM(i, :)';
            Cinv_i = precision_lnorm(:, :, i);
            summands(:, :, i) = inv(Cinv_i+Dinv) + (beta_i-b') * (beta_i-b')';
        end
        D = mean(summands, 3);                                              % update D
        
        if max(eucl_rel(b, b_old), eucl_rel(D, D_old)) < 1e-3
            break                                                           % break at convergence tolerance
        end
    end
end


function precision = invert_uncertainties(varbeta_lnorm)                % Inverse covariance matrix estimates
    precision = zeros(size(varbeta_lnorm));
    for i = 1:size(varbeta_lnorm, 3)
        covariance = varbeta_lnorm(:, :, i);

        fixed = diag(covariance) < max(diag(covariance)) / 1e6;             % manual correction for near-zero uncertainty
        precision(~fixed, ~fixed, i) = tryinv(covariance(~fixed, ~fixed));
        precision(fixed, fixed, i) = diag(Inf(1, sum(fixed)));
    end
end


function d = eucl_rel(x, y, ~)                                          % Relative Euclidean distance
    if nargin < 3
        if iscell(x) && iscell(y)
            x = cell2mat(x);
            y = cell2mat(y);
        end
        x = flatten(x);
        y = flatten(y);
    end
    
    d = norm(x - y) / norm(x);
end


function Ainv = tryinv(A)
    try
        Ainv = inv(A);
    catch ME
        if strcmp(ME.identifier, 'MATLAB:nearlySingularMatrix')
            Ainv = pinv(A);
        end
    end
end