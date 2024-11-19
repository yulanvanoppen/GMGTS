function d = hedist(m1, S1, m2, S2)
    m1 = flatten(m1);
    m2 = flatten(m2);
    
    d = sqrt(1 - 2^(length(m1)/2) * det(S1*S2)^.25 * det(S1+S2)^-.5 ...
                * exp(-1/4 * (m1-m2)' * inv(S1 + S2) * (m1-m2)));
end