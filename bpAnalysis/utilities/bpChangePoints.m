function cp = bpChangePoints(Data, dim, nShuff, direction)
    
if nargin < 2
    dim = 1;
end

if nargin < 3
    nShuff = 1000;
end

if nargin < 4
    direction = 'both'; % choices; 'up', 'down', or 'both'[default] 
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
    
    null = zeros(nShuff, 1);
    
    %% see comments below further down on rationale for this code block (better for speed, but hacky)
    switch direction
        case 'both'
            [val, ix] = changepoint_both(data(valid), nochange);
            for shuffCounter = 1:nShuff % null distribution
                shuffIx = randperm(np);
                [nullVal, ~] = changepoint_both(data(valid(shuffIx)), nochange);
                null(shuffCounter) = nullVal;
            end
        case 'up'
            [val, ix] = changepoint_up(data(valid), nochange);
            for shuffCounter = 1:nShuff % null distribution
                shuffIx = randperm(np);
                [nullVal, ~] = changepoint_up(data(valid(shuffIx)), nochange);
                null(shuffCounter) = nullVal;
            end            
        case 'down'
            [val, ix] = changepoint_down(data(valid), nochange);
            for shuffCounter = 1:nShuff % null distribution
                shuffIx = randperm(np);
                [nullVal, ~] = changepoint_down(data(valid(shuffIx)), nochange);
                null(shuffCounter) = nullVal;
            end            
    end
    %%
    


    difference = abs(nochange - cumsum(data(valid))); % redundant to subfunction
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


% Hack for speed- I created 3 versions of subfunction to avoid use of
% switch condition (a string) in subfunction- very slow

% better way to do this would be to delegate null distribution generation
% to subfunction (one call vs nShuffles calls to subfunction)

function [val, ix] = changepoint_both(data, nochange)
test = cumsum(data);
difference = nochange - test;
difference = abs(difference);
[~, ix] = max(difference);
val = abs(difference(ix));


function [val, ix] = changepoint_up(data, nochange)
test = cumsum(data);
difference = nochange - test;
[~, ix] = max(difference);
val = abs(difference(ix));

function [val, ix] = changepoint_down(data, nochange)
test = cumsum(data);
difference = nochange - test;
[~, ix] = min(difference);
val = abs(difference(ix));    
    
% THIS IS SLOWER:    
% function [val, ix] = changepoint_both(data, nochange, direction)
% test = cumsum(data);
% difference = nochange - test;
% switch direction
%     case 'both'
%         difference = abs(difference);
%         [~, ix] = max(difference);
%     case 'up'
%         [~, ix] = max(difference);
%     case 'down'
%         [~, ix] = min(difference);
% end
% val = abs(difference(ix));


            
        
    