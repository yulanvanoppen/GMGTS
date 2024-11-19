function [estimates, estimator] = GMGTS(model_file, data, varargin)
%GMGTS Infer random effects distribution using GMGTS algorithm
%   Detailed explanation goes here

    system = ODEIQM(model_file, varargin{:}); 
    estimator = Estimator(system, data, varargin{:});
    
    outputs = cell(1, length(estimator.method));
    [outputs{:}] = estimator.estimate();
    
    if length(estimator.method) == 1
        estimates = outputs{1};
    else
        estimates = [outputs{1} outputs{2}];
    end
    
    plot(estimator)
end

