function d = eucl_rel(x, y, ~)
    if nargin < 3
        x = reshape(x, [], 1);
        y = reshape(y, [], 1);
    end
    
    d = norm(x - y) / norm(x);
end