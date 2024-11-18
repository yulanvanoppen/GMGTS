function varargout = mat2csl(mat)
    mat = mat2cell(mat, size(mat, 1), ones(1, size(mat, 2)));
    varargout = mat{:};
end

