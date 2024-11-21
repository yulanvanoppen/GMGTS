function [out, estimator] = GMGTS(model_file, data, varargin)
%GMGTS infers random effects distributions of ODE-based ME models.
%   GMGTS attempts to recover the distribution parameters b,D of a
%   (linear in parameters) ODE-based mixed-effects model
%       dx_i(t)/dt = g(t, x_i(t)) beta_i + h(t, x_i(t))
%       beta_i ~ (log)Normal(b, D)                          (i = 1, ..., N)
%   from measurements (t_j, y_ij), where y_ij are vectors of observed
%   components of x_i(t_j) perturbed by measurement noise.
%
%   GMGTS builds on the GTS framework, using gradient matching to obtain
%   cell-specific estimates after smoothing the measurements.
%
%   out = GMGTS(model_file, data, ...) infers random effect
%   distributions of the system specified in the IQM model_file
%   (instructions for setting up the model file are given below) from
%   measurements given in data. Here data is either a TxLxN-dimensional
%   array with measurements at T time points for L observables and N cells
%   (individuals), or a 1x1 struct with its measurements stored in a field
%   named y or traces. The estimates are returned as a struct containing
%   the inferred random effect mean b and covariance matrix D, individual
%   estimates beta, and predicted states. Additional arguments are passed
%   to the System and Estimator constructors, as well as Estimator.plot() 
%   method, see the details below.
% 
%   out = GMGTS(model_file, data, t, ...) assumes which the
%   measurements were taken at time points t. If data is a struct, t is
%   ignored and assumed to be a field of data. If omitted, the default
%   t = 0:size(data, 1)-1 is used.
% 
%   out = GMGTS(model_file, data, t, observed, ...) specifies the
%   indices of the observables with respect to the system determined by
%   model_file through observed. If data is a struct, observed is ignored
%   and assumed to be a field of data. If omitted, it is assumed that
%   observed = 1:size(data, 2).
% 
%   out = GMGTS(model_file, data, t, observed, init, ...) integrates
%   the ODE system from the initial values given in init to make state
%   predictions. If data is a struct, init is ignored and assumed to be a
%   field of data. If omitted, the initial values are all assumed to equal
%   1e-8.
%
%   out = GMGTS(model_file, data, t, observed, init, Plot=false, ...)
%   disables plots with parameter estimates, the inferred random effects 
%   distribution, model predictions, and any smoothed measurements (enabled
%   by default).
% 
%   [out, estimator] = GMGTS(model_file, data, ...) also returns the
%   instantiated Estimator object.
%
%   See also SYSTEM, ESTIMATOR.

    default_Plot = true;                                                    % add argument whether or not to plot
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'Plot', default_Plot, @islogical);
    parser.parse(varargin{:})

    system = ODEIQM(model_file, varargin{:});                               % create (ODE) System object
    estimator = Estimator(system, data, varargin{:});                       % create Estimator object from system and data
    
    outputs = cell(1, length(estimator.method));
    [outputs{:}] = estimator.estimate();                                    % collect estimates
    
    if length(estimator.method) == 1                                        % concatenate estimates if both GMGTS and GTS used
        out = outputs{1};
    else
        out = [outputs{1} outputs{2}];
    end
    
    if parser.Results.Plot, plot(estimator, varargin{:}), end               % plot if required
end

