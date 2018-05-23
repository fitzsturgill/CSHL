function out = stripNaNs(in)
    % strips NaNs from vectors
    
    out = in(~isnan(in));