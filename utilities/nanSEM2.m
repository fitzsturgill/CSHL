function s = nanSEM2(A, dim)
    % pointwise computation of SEM, takes matrix as input
    if nargin < 2
        dim = 1;
    end
    
    s = nanstd(A, 0, dim) ./ sqrt(sum(isfinite(A), dim));
    