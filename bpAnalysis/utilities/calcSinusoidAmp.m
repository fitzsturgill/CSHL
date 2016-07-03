function [amp minValue maxValue] = calcSinusoidAmp(modData, nBins, show)
    % estimate amplitude of a (somewhat) noisy sinusoid based commonly occuring minimum and
    % maximum values
    if nargin < 2
        nBins = 200; % number of bins across which to search for common minimum and maximum values of each cycle
    end
    
    if nargin < 3
        show = 0; % to validate based upon data and choice of nBins
    end
    
    % generate a histogram of the values in the periodic data (the mins and
    % maxes are more common due to their shallow slope)
    
    [N, edges] = histcounts(modData, nBins);
    centers = (edges(1:end-1) + edges(2:end)) /2;
    
    % split it in half (this depends on the U-shaped (paul says Beta)
    % distribution of y values
    
    splitPoint = floor(nBins/2);
    
    [~, minPoint] = max(N(1:splitPoint)); % most frequently occuring bin in lower half of data
    minValue = centers(minPoint);
    
    [~, maxPoint] = max(N(splitPoint + 1 : end));
    maxPoint = maxPoint + splitPoint;
    maxValue = centers(maxPoint);
    
    amp = maxValue - minValue;
    
    if show
        h = ensureFigure('calcSinusoidAmp', 1);
        plot(centers, N, 'k'); hold on
        plot([minValue minValue], [0 N(minPoint)]);
        plot([maxValue maxValue], [0 N(maxPoint)]);
        ylabel('Histogram Frequency');
        xlabel('Y axis values (bins centers)');
    end

    
    
    
    
    
    