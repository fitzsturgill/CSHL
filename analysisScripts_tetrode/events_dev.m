            
%%
pname = 'C:\FitzData\Cellbase_dev\CD_3\170821a\';
% pname = 'Z:\FitzTetrode\';
fname = 'Events.nev';new_fname='EVENTS';
param2 = [1 1 1 1 1];
param3 = 1;
param4 = 1;
param5 = [];
[Events_TimeStamps, Events_EventIDs, Events_Nttls, Events_Extras, Events_EventStrings, Events_NlxHeader] = Nlx2MatEV([pname fname],param2,param3,param4,param5);

%%
stimN = 128;
stimIdx = Events_Nttls == stimN;
stim_IDs = Events_EventStrings(stimIdx);

%u
%% 

fitz = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0080).';
expr = '(?<=port.).(?=.value)';
portS = regexp(fitz, expr, 'match');

% event 128:
% parsed_ttls(2,:)
% 
% ans =
% 
%      1     0     0     0     0     0     0     0

%% find the most common event (should be the laser TTL event)

uniqueEvents = unique(Events_Nttls);
eventCounts = zeros(size(uniqueEvents));

for counter = 1:length(uniqueEvents)
    eventCounts(counter) = sum(Events_Nttls == uniqueEvents(counter));
end
[counts_sorted, counts_sortIndices] = sort(eventCounts);
events_sorted = uniqueEvents(counts_sortIndices);

%% extract the binary representations

events_bin = dec2bin(events_sorted);
events_binTrunc = events_bin(:, 1:4);
events_trunc = bin2dec(events_binTrunc);

%% ITIs of event 128

laser_idx = Events_Nttls == 128;
laser_times = Events_TimeStamps(laser_idx);
laser_itis = diff(laser_times);

ensureFigure('laser_times', 1);
stem(laser_times, ones(size(laser_times)));






