function dPrime = dPrime_SNR(smallerMean,largerMean)
% difference in means divided by square root of mean variance
    dPrime = nanmean(largerMean) - nanmean(smallerMean)...
        / (sqrt((nanstd(smallerMean)^2 + nanstd(largerMean)^2)/2)); 