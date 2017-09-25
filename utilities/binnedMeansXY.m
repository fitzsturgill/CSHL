function [meansY, semY, binCenters] = binnedMeansXY(x, y, nbins)


    bins = linspace(min(min(x)), max(max(x)), nbins + 1);
    binCenters = bins(1:end-1) + bins(2) - bins(1); 
    meansY = zeros(nbins, 1);
    semY = zeros(nbins, 1);
    for counter = 1:nbins
        inThisBin = x >= bins(counter) & x < bins(counter + 1);
        meansY(counter) = nanmean(y(inThisBin));
        semY(counter) = nanSEM(y(inThisBin));
    end
    
    