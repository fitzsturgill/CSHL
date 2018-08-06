% test_script for 
% loadcb; cellid=CELLIDLIST{1}; tagging(cellid); taggedpropTO(cellid); taggedsummaryTO(cellid)

% figure; viewcell2b(thiscell,'TriggerName','Us_start','SortEvent','trialNumber','eventtype','behav','ShowEvents',{'Us_start'},...
%    'Partitions','#trialType: {1 4 7}','window',[-7 4], 'dt', 0.01, 'sigma', 0.02, 'PSTHstd', 'on', 'isadaptive', true);

window = [-7 4];
dt = 0.01;
sigma = 0.02;
winMargin = sigma * 3;
rat = 'CP9';
session = '180731a';
thesecells = findcell('rat', rat, 'session', session);

loadcb; % % loads ANALYSES, CELLIDLIST, PREFERENCES, and TheMatrix
pref = getcbpref; % torben's new way

cellid = thesecells{4};
TE = loadcb(cellid, 'TrialEvents');
SP = loadcb(cellid,'EVENTSPIKES');


%% add number of outcome Licks occuring between 0.1 and 1s following outcome
TE.usLicks = countEventFromTE(TE, 'Port1In', [0.1 1], TE.Us);
TE.usLicks_rate = TE.usLicks.rate;
TE.csLicks = countEventFromTE(TE, 'Port1In', [0 3], TE.Cue);
TE.csLicks_rate = TE.csLicks.rate;
TE.rewardReceived = TE.usLicks.count >= 2; % get rid of trials where mouse fails to lick for reward
save(fullfile(pref.datapath, rat, session, pref.TrialEvents_fname), 'TE'); 


% make rasters and PSTHs
time = window(1) - winMargin:dt:window(2) + winMargin;
binraster = stimes2binraster(stimes,time,dt);
%%
ensureFigure('test1', 1); viewcell2b(cellid,'TriggerName','Us_start','SortEvent','usLicks_rate','eventtype','behav','ShowEvents',{'Us_start'},...
   'Partitions','#trialType: {1 4 7} & rewardReceived','window',[-7 4], 'dt', 0.01, 'sigma', 0.02, 'PSTHstd', 'on', 'isadaptive', true);
%%
ensureFigure('test2', 1); viewcell2b(cellid,'TriggerName','Us_start','SortEvent','csLicks_rate','eventtype','behav','ShowEvents',{'Us_start'},...
   'Partitions','#trialType: {1 4 7}','window',[-7 4], 'dt', 0.02, 'sigma', 0.02, 'PSTHstd', 'on', 'isadaptive', true);
%%
ensureFigure('test3', 1); viewcell2b(cellid,'TriggerName','Us_start','SortEvent','csLicks_rate','eventtype','behav','ShowEvents',{'Us_start'},...
   'Partitions','all','window',[-7 4], 'dt', 0.02, 'sigma', 0.02, 'PSTHstd', 'on', 'isadaptive', true);


