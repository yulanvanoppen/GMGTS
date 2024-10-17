function d = eucl_rel(x, y, ~)
    if nargin < 3
        if iscell(x) && iscell(y)
            x = cell2mat(x);
            y = cell2mat(y);
        end
        x = reshape(x, [], 1);
        y = reshape(y, [], 1);
    end
    
    d = norm(x - y) / norm(x);
end