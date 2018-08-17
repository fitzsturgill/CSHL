% load(fullfile('F:\Cellbase\CP6\180703a', 'Events.mat'));   % load Neuralynx events
% TrialStart_nlx = getBehaviorStartTimes(Events_Nttls, Events_EventStrings, Events_TimeStamps);
% 
nChunks = size(sampleCheck.beforeSegments, 1);
answers = zeros(nChunks, 3);
range = 1:100;
for counter = 1:nChunks
    answers(counter, 1) = strfind(sampleCheck.samples', sampleCheck.beforeSegments(counter, range));
    answers(counter, 2) = strfind(sampleCheck.samples', sampleCheck.afterSegments(counter, range));
    try
        answers(counter, 3) = strfind(sampleCheck.samples', [sampleCheck.afterSegments(counter, :) sampleCheck.beforeSegments(min(counter + 1, nChunks), :)]);
    catch
        answers(counter,3) = 0;
    end
end

% 
% %%
% starts = find(Events_Nttls == 128);
% nexts = starts + 1;
% 
% 
% 
% %%
% % function TrialStart_nlx = trytomatch(ts_bpod,TrialStart_nlx)
% function TrialStart_nlx = trytomatch(ts,TrialStart_nlx)
% 
% % Try to match time series by shifting
% len = length(TrialStart_nlx) - 15;
% minx = nan(1,len);
% for k = 1:len
%     minx(k) = max(diff(ts(1:15)-TrialStart_nlx(k:k+14)));  % calculate difference in the function of shift
% end
% mn = min(abs(minx));
% minx2 = find(abs(minx)==mn);
% minx2 = minx2(1);   % find minimal difference = optimal shift
% TrialStart_nlx = TrialStart_nlx(minx2:min(minx2+length(ts)-1,length(TrialStart_nlx)));