            
%%
pname = 'Z:\FitzTetrode\2017-07-13_16-20-52\';
% pname = 'Z:\FitzTetrode\';
fname = 'Events.nev';new_fname='EVENTS';
param2 = [1 1 1 1 1];
param3 = 1;
param4 = 1;
param5 = [];
[Events_TimeStamps, Events_EventIDs, Events_Nttls, Events_Extras, Events_EventStrings, Events_NlxHeader] = Nlx2MatEV([pname fname],param2,param3,param4,param5);

%% find the most common event (should be the laser TTL event)

uniqueEvents = unique(Events_Nttls);
eventCounts = zeros(size(uniqueEvents));

for counter = 1:length(uniqueEvents)
    eventCounts(counter) = sum(Events_Nttls == uniqueEvents(counter));
end
[counts_sorted, counts_sortIndices] = sort(eventCounts);
events_sorted = uniqueEvents(counts_sortIndices);


