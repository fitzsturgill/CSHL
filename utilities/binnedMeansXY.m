function [meansY, semY, binCenters] = binnedMeansXY(x, y, binSpec)
    % bins- can be supplied as a scalar to specify the total number of bins
    % or as a vector corresponding to the bin edges, can be unequal in size

    if isscalar(binSpec)
        bins = linspace(min(min(x(isfinite(x)))), max(max(x(isfinite(x)))), binSpec + 1);
    else
        bins = binSpec;
    end
    
    binCenters = (bins(1:end-1) + bins(2:end))/2; 
    meansY = zeros(length(binCenters), 1);
    semY = zeros(length(binCenters), 1);
    for counter = 1:length(binCenters)
        inThisBin = x >= bins(counter) & x < bins(counter + 1);
        meansY(counter) = nanmean(y(inThisBin));
        semY(counter) = nanSEM(y(inThisBin));
    end
    
    