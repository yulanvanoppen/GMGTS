function Ainv = tryinv(A)
    try
        Ainv = inv(A);
    catch ME
        if strcmp(ME.identifier, 'MATLAB:nearlySingularMatrix')
            Ainv = pinv(A);
        end
    end
end