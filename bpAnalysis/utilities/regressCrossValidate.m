
function Rsq = regressCrossValidate(y, X, nFolds)
    if nargin < 3
        nFolds = 10;
    end
    n = length(y);
    shortIx = round(n/10); % end of the decile corresponding to left out 10% of data
    predicted = zeros(shortIx, nFolds); % accumulated predictions applied to held out 10% of data
    target = zeros(shortIx, nFolds); % matching target values corresponding to held out 10% of data
    for iFolds = 1:nFolds
        idx = randperm(length(y));
        idx = idx(1:shortIx);
        cIdx = false(n, 1);
        cIdx(idx) = true; % create logical index to access random 10% of data
        
        theseWeights = regress(y(~cIdx),X(~cIdx, :)); %full model
        predicted(:,iFolds) = X(cIdx, :) * theseWeights;
        target(:, iFolds) = y(cIdx);
    end
    
    Rsq = corr2(target(:), predicted(:));
    
        
        
        