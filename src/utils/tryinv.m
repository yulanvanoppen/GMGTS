function Ainv = tryinv(A)
    try
        Ainv = inv(A);
    catch ME                                                                % fall back on the pseudoinverse for
        if strcmp(ME.identifier, 'MATLAB:nearlySingularMatrix')             % (nearly) singular matrices when the
            Ainv = pinv(A);                                                 % corresponding warning is set to 'error'
        else
            rethrow(ME)
        end
    end
end