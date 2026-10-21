%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create latex table from matlab table
% 
% @author Ing. Ales Jandera
% @email ales.jandera@tuke.sk
%
% MIT license
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = table2latex(varargin)
    % define initial variables
    numArgs = nargin;
    cr = 13;
    tab = 9;
    format = '%0.2f';
    N = length(varargin{1});
    M = numArgs;
    align = 'l';
    
    % validation of input arguments
    if numArgs < 1
        error('Not enough input arguments. Please provide at least one.');
    end
    
    % parse format if inserted
    for m = 1:numArgs
       name = inputname(m);
       % parse format parameters
       if name == "format"
           format = varargin{m};
           M = M - 1;
       elseif name == "align"
           align = varargin{m};
           M = M - 1;
       end
    end
    % begin the latext table
    str = ['\begin{table}', cr, tab, '\centering', cr, tab, tab,...
        '\begin{tabular}{', repmat(align,1,M),'}', cr, tab, tab, tab];
    % Work through header elements
    for n = 1:N-1
       name = inputname(n);
       if ischar(name)
           str = [str, '\textbf{', name,'} & '];
       else
           str = [str, '\textbf{',num2str(name,format), '} & '];
       end
    end
    
    str = [str, '\textbf{', num2str(inputname(n+1), format), '}\\',...
        cr,tab,tab,tab];
    
    % Work through table elements
    for n = 1:N
        for m = 1:M-1
           value = varargin{m}(n);
           if ischar(value)
               str = [str, '$', value, '$ & '];
           else
               str = [str, '$', num2str(value,format), '$ & '];
           end
        end
    
        str = [str, '$', num2str(varargin{m+1}(n),format), '$ \\',...
            cr,tab,tab,tab];
    end
    
    % remove unusable end of the string
    str = str(1:end-6);
    
    % add close part of the latex table
    str = [str, cr, tab, tab '\end{tabular}', cr, '\end{table}'];
    
    if nargout == 0
        % write to clipboard
        clipboard('copy', str);
    else
        % write to output argument
        varargout{1} = str;
    end
end
