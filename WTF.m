
% find trials to 50%
fraction = 0.5;
compFields = {'csPlus', 'csMinus';...  % first row is after reversal
              'csMinus', 'csPlus'...   % second row is before reversal
              };
fitFields = {'licks_cs', 'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2'};
outputFields = {'bl', 'max', 'tt50'};
tt50 = struct();

nReversals = size(AR.csPlus.globalTrialNumber.after, 1);

for counter = 1:length(compFields)
    for counter2 = 1:length(fitFields)
        for counter3 = 1:length(outputFields)            
            tt50.(compFields{counter}).(fitFields{counter2}).(outputFields{counter3}) = NaN(nReversals, 1);
        end
    end
end
tt50.fraction = fraction;
tt50.baselineTrials = baselineTrials;


for compCounter = 1:size(compFields, 2)    
    for fieldCounter = 1:length(fitFields)
        bl = nanmean(AR.(compFields{2, compCounter}).(fitFields{fieldCounter}).before(:, end - baselineTrials + 1:end), 2);
        tt50.(compFields{1, compCounter}).(fitFields{fieldCounter}).bl = bl;
        tt50Data = AR.(compFields{1, compCounter}).(fitFields{fieldCounter}).after;
        for counter = 1:size(tt50Data, 1)        
            thisData = tt50Data(counter, :);
            if ~any(isfinite(thisData))
                continue
            end
            switch compFields{1, compCounter}
                case 'csPlus'
                    top = percentile(thisData, 0.9);
                    thresh = (bl(counter) + (top - bl(counter)) * fraction);
                    latency = find(thisData > thresh, 1);
                    if isempty(latency)
                        continue
                    end
                    tt50.(compFields{1, compCounter}).(fitFields{fieldCounter}).tt50(counter) = latency;
                    tt50.(compFields{1, compCounter}).(fitFields{fieldCounter}).thresh(counter) = thresh;
                case 'csMinus'
                    bottom = percentile(thisData, 0.1);
                    thresh = (bl(counter) - (bl(counter) - bottom) * fraction);
                    latency = find(thisData < thresh, 1);
                    if isempty(latency)
                        continue
                    end
                    tt50.(compFields{1, compCounter}).(fitFields{fieldCounter}).tt50(counter) = latency;                    
                    tt50.(compFields{1, compCounter}).(fitFields{fieldCounter}).thresh(counter) = thresh;
            end            
        end
    end
end


%%

% dT = 0.05;
% 
% sr = 40;
% 
% duration = 60;
% 
% nTrials = 50;
% 
% stimes = cell(nTrials, 1);
% for counter = 1:nTrials
%     stimes{counter} = makeSpikes(dT, sr, duration, 1);
% end
% fitz= 
% binraster = stimes2binraster(stimes, 0:dT:duration, dT, repmat([0 60], nTrials, 1));
% 
% 
% 
% % binraster = stimes2binraster(spikes
% 
% 
% function spikes = makeSpikes(timeStepS, spikesPerS, durationS, numTrains)
% 
% if (nargin < 4)
%     numTrains = 1;
% end
% times = [0:timeStepS:durationS];
% spikes = zeros(numTrains, length(times));
% for train = 1:numTrains
%     vt = rand(size(times));
%     spikes(train, :) = (spikesPerS*timeStepS) > vt;
% end
% end