
cp.csPlus.licks_cs = bpChangePoints([AR.csMinus.licks_cs.before(:, end - baselineTrials + 1:end) AR.csPlus.licks_cs.after(:, 1:end)], 2, 1000, 'up');


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