%% optogenetic tagging graphs for research statement
cellid = findcell('rat', 'CD4', 'session', '170909a');
TE = loadcb(cellid,'StimEvents');
SP = loadcb(cellid,'STIMSPIKES');

% View light-triggered raster and PSTH
TrigEvent = 'Pulse';
SEvent = 'Pulse';
win = [-0.05 0.05];
parts = 'all';
%     parts = '#BurstNPulse';
dt = 0.001;
sigma = 0.001;
PSTHstd = 'on';
ShEvent = {'Pulse'};
ShEvColors = hsv(length(ShEvent{1}));
ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);

% viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',ShEvColors,...
%     'FigureNum',H,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
%     'EventMarkerWidth',0,'PlotZeroLine','off')    

%%
trigger_pos = findcellstr(SP.events(:,1),'Pulse');
TriggerEvent = SP.events{trigger_pos,2};
stimes  = SP.event_stimes{trigger_pos};

window_margin = SP.events{trigger_pos,4};
ev_windows = SP.event_windows{trigger_pos};
time = win(1):dt:win(2);
binraster = stimes2binraster(stimes,time,dt);
trials = length(stimes);
% valid_trials = find(~isnan(TE.(TriggerEvent)));
valid_trials = trials; % try including trials without a spike to give sense of reliability in raster

[psth, spsth, spsth_se] = binraster2psth(binraster,dt,sigma,COMPTRIALS,valid_trials);
% trials = length(