function setup_GMGTS(save_path)
%SETUP_GMGTS adds the GMGTS code, models, and data folders to the MATLAB path.
%
%   SETUP_GMGTS() adds the source code (src/), model files (models/), data
%   (data/), and experiment utilities (experiments/utils/) to the MATLAB
%   path for the current session and checks for the required toolboxes.
%   Run this once per MATLAB session from any folder before using GMGTS,
%   running the examples, or reproducing the results in the paper.
%
%   SETUP_GMGTS(true) additionally saves the path for future sessions.
%
%   See also GMGTS, SYSTEM, GENERATOR, ESTIMATOR.

    if nargin < 1, save_path = false; end
    root = fileparts(mfilename('fullpath'));

    addpath(fullfile(root, 'src'), fullfile(root, 'src', 'utils'), ...
            fullfile(root, 'models'), genpath(fullfile(root, 'data')), ...
            fullfile(root, 'experiments', 'utils'), fullfile(root, 'examples'));
    fprintf('GMGTS folders added to the MATLAB path (root: %s).\n', root);

    required = {'Symbolic Math Toolbox', 'Optimization Toolbox', ...        % used by System, Optimization, and
                'Statistics and Machine Learning Toolbox', ...              % random effect sampling/plots
                'Model Predictive Control Toolbox'};                        % (mpcActiveSetSolver)
    toolboxes = ver;
    installed = {toolboxes.Name};
    missing = required(~ismember(required, installed));
    if ~isempty(missing)
        warning('GMGTS:missingToolbox', 'Missing MATLAB toolbox(es): %s.', strjoin(missing, ', '));
    end
    if exist('installIQMtools', 'file') ~= 2
        warning('GMGTS:missingIQM', ['IQM Tools (installIQMtools) not found on the path; ' ...
                'constructing System objects from model files will fail. ' ...
                'See README.md for installation instructions.']);
    end

    if save_path, savepath; end
end
