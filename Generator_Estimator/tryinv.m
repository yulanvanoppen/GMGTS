function Ainv = tryinv(A)
    try
        Ainv = inv(A);
    catch ME
        if strcmp(ME.identifier, 'MATLAB:nearlySingularMatrix')
            Ainv = pinv(A);
            disp(111111111111111111111111111111111111111111111111111111111111)
        end
    end
end