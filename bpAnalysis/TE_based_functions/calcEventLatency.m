
function lags = calcEventLatency(TE, eventField, zeroTimes, endTimes)
%  inputs: TE, eventField = TE field that contains cell arrays of event
%  times (e.g. Port1In), zeroTimes- vector or cell array of zero times, if
%  cell array, takes the first element in it as the zero time for that
%  trial..., endTimes- vector or cell array of times corresponding to the
%  end of the time window in which you are looking for the first event and
%  its corresponding latency.  If event occurs after endtimes

    if nargin < 4
        endTimes = [];
    end

    nTrials = length(TE.(eventField));
    lags = nan(nTrials, 1);
    for counter = 1:nTrials
        theseEvents = TE.(eventField){counter};
        if iscell(zeroTimes)
            thisZero = zeroTimes{counter}(1);
        else
            thisZero = zeroTimes(counter);
        end
        if iscell(endTimes)
            thisEnd = endTimes{counter}(1);
        elseif isempty(endTimes)
            thisEnd = Inf;
        else
            thisEnd = endTimes(counter);
        end        
        first = find((theseEvents > thisZero) & (theseEvents <= thisEnd), 1);
        if ~isempty(first)
            lags(counter) = theseEvents(first) - thisZero;
        end
    end