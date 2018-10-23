function cp = bpChangePoints(Data, dim, nShuff)
    
if nargin < 2
    dim = 1;
end

if nargin < 3
    nShuff = 1000;
end

assert(ndims(Data) <= 2);
if dim == 2
    Data = Data';
end

ncp = size(Data, 2);
cp = struct(...
    'index', NaN(ncp, 1),...
    'p', NaN(ncp, 1),...
    'logit', NaN(ncp, 1)...
    );
cp.cumsum = cell(ncp, 1);
cp.logitAll = cell(ncp, 1);
cp.distance = cell(ncp, 1);

for counter = 1:ncp
    data = Data(:,counter);
    valid = find(~isnan(data));
    np = length(valid);
    total = sum(data(valid));
    nochange = total/np * (1:np); nochange = nochange(:); % diagonal line from orgin
    [val, ix] = changepoint(data(valid), nochange);
    null = zeros(nShuff, 1);
    for shuffCounter = 1:nShuff % null distribution
        shuffIx = randperm(np);
        [nullVal, ~] = changepoint(data(valid(shuffIx)), nochange);
        null(shuffCounter) = nullVal;
    end
    difference = nochange - cumsum(data(valid)); % redundant to subfunction
    p = iprctile(null, val);
    cp.index(counter) = valid(ix);
    cp.p(counter) = p;
    cp.logit(counter) = log(p/(1-p));
    cp.cumsum{counter} = cumsum(data(valid));
    pAll = zeros(np, 1);
    for point = 1:np
        pAll(point) = iprctile(null, difference(point));
    end
    cp.logitAll{counter} = log(pAll ./ (1 - pAll));
    cp.distance{counter} = difference / total;
end




    
    
    
    
function [val, ix] = changepoint(data, nochange)
test = cumsum(data);
difference = nochange - test;
[~, ix] = max(difference);
val = difference(ix);


            
        
    