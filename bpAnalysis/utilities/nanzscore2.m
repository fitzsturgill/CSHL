function [z,mu,sigma] = nanzscore2(x)

    x2 = reshape(x, numel(x), 1);
    mu = nanmean(x2);
    sigma = nanstd(x2);
    z = (x - mu) / sigma;