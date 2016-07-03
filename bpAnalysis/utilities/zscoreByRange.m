function [Z mu sigma] = zscoreByRange(X, i1, i2)
    % standardize vector by a portion of it
    [~, mu, sigma] = zscore(X(i1:i2));
    Z = X - mu;
    Z = Z / sigma;
    
    