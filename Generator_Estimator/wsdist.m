function d = wsdist(m1, S1, m2, S2)
    m1 = reshape(m1, [], 1);
    m2 = reshape(m2, [], 1);
    
    s2 = sqrtm(S2);
    
    WS = sqrt(sum((m1-m2).^2) + trace(S1 + S2 - 2*sqrtm(s2 * S1 * s2')));
    
    d = WS;
end