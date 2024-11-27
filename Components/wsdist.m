function d = wsdist(m1, S1, m2, S2, normalized)
    if nargin < 5, normalized = true; end
    m1 = flatten(m1);
    m2 = flatten(m2);
    
    if normalized
        S1 = S1 ./ m1 ./ m1';
        S2 = S2 ./ m1 ./ m1';
        m2 = m2 ./ m1;
        m1 = m1 ./ m1;
    end
    s2 = sqrtm(S2);
    
    d = sqrt(sum((m1-m2).^2) + trace(S1 + S2 - 2*sqrtm(s2 * S1 * s2')));
end