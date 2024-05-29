function d = hedist(m1, S1, m2, S2)
    m1 = reshape(m1, [], 1);
    m2 = reshape(m2, [], 1);
    
    HE = sqrt(1 - 2^(length(m1)/2) * det(S1*S2)^.25 * det(S1+S2)^-.5 * exp(-1/4 * (m1-m2)' * inv(S1 + S2) * (m1-m2)));
    
    d = HE;
end