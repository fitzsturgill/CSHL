function [R, shiftR, rawR, lags] = correctedXCorr(x, y, maxlag, dim)

% computes the average, all-ways shift-predictor subtracted (corrected)
% cross-correlogram between matrices x, y
%     if dim == 1 data lies in columns, if dim == 2, rows
    if nargin < 4
        dim = 1;
    end
    
    if ndims(x) > 2
        error('');
    end
    
    if ~all(size(x) == size(y))
        error('');
    end
    
    if dim == 2
        x = x';
        y = y';
    end
    
    if nargin < 3 || isempty(maxlag)
        maxlag = size(x,1) - 1;
    end
    
    x = nanzscore(x, 0);
    y = nanzscore(y, 0);
    
    [rawR, lags] = avgXCorr(x, y, maxlag); % uncorrected

    %     shift predictor
    shiftR = zeros(maxlag*2+1, size(x,2)-1);
    for counter = 1:size(x,2)-1
        sy = circshift(y, [0 counter]);
        [sy, ~] = avgXCorr(x, sy, maxlag);
        shiftR(:,counter) = sy;
    end
    shiftR = nanmean(shiftR, 2);
    R = rawR - shiftR;

end



function [R, lags] = avgXCorr(x, y, maxlag)
    allR = zeros(maxlag*2+1, size(x,2));
    lags = -maxlag:maxlag;

    for column = 1:size(x, 2)
        [theseR, lags] = xcorr(x(:,column), y(:,column), maxlag);
        allR(:,column) = theseR;
    end
    R = nanmean(allR, 2);
end
    
