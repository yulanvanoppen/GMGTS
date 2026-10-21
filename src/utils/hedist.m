function d = hedist(m1, S1, m2, S2)
%HEDIST Hellinger distance between two multivariate normal distributions.
%
%   d = HEDIST(m1, S1, m2, S2) returns the Hellinger distance between
%   N(m1, S1) and N(m2, S2), defined through
%       d^2 = 1 - det(S1)^(1/4) det(S2)^(1/4) / det(S)^(1/2)
%                 * exp(-1/8 (m1-m2)' S^(-1) (m1-m2)),   S = (S1+S2)/2,
%   so that 0 <= d <= 1. The distance is invariant under invertible affine
%   transformations, so no normalization (cf. WSDIST) is applied.
%
%   See also WSDIST.

    m1 = reshape(m1, [], 1);
    m2 = reshape(m2, [], 1);
    S1 = (S1 + S1') / 2;
    S2 = (S2 + S2') / 2;
    S = (S1 + S2) / 2;

    logdet = @(A) sum(log(max(eig(A), realmin)));                          % stable log-determinants
    log_bc = logdet(S1)/4 + logdet(S2)/4 - logdet(S)/2 ...                  % log Bhattacharyya coefficient
             - (m1 - m2)' * (S \ (m1 - m2)) / 8;

    d = sqrt(max(0, 1 - exp(log_bc)));
end
