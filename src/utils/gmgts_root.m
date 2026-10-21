function root = gmgts_root()
%GMGTS_ROOT returns the root folder of the GMGTS repository.
%
%   Used by the scripts in examples/ and experiments/ to locate data files
%   and to store results independently of the current folder.

    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
