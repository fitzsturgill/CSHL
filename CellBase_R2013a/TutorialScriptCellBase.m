% Note that I've placed P35, P36 in CSHL repo so that I have example data
% for the future that runs with the TutorialScript.  Don't initialize
% cellbase in the repo



%%%     EXAMPLE SCRIPT FOR CELLBASE     %%%
%%%                                     %%%
%%%         AK 11/2006  PM 07/2017                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the dabase -- needs to be run once -- fast.
%
initcb

%%
% Using spike times & event times prealigns and saves the data according
% several criteria  -- run once, will take ~15 minutes.
%

% Define events of interest
[events,epochs] = defineEventsEpochs_Dual2AFCv3;

%%
prealignSpikesPM2('all','FUNdefineEventsEpochs',@defineEventsEpochs_Dual2AFCv3,'filetype','event',...
    'events',[],'epochs',[],'writing_behavior','overwrite','ifsave',1);

%%
allcells = listtag('cells');       % list all cellids
%%
% Let's add an analysis 
% Runs an analysis on all units using the function 'meanrate'
addanalysis(@meanratePM,'property_names',{'MeanRate'});

%Lists all the known properties, results of analyses
listtag('properties')   % or it's enough to type: listtag('prop')
%%
% List the value of the property 'MeanRate' for all the units
rates = getvalue('MeanRate');

% Access data from a cell list

cellid=allcells{1};

TE = loadcb(cellid,'TrialEvents');  % Behavioral data
SP = loadcb(cellid,'EVENTSPIKES'); % Spiking data

% Lets plot a psth


[TRIALS, TAGS, RTAGS, NUM_TAGS, PartNum] = partition_trials(TE,'#SideReward');

% Trigger on start of anticipation period
TriggerEvent = 'ResponseStart';

% SortingEvent='ResponseEnd';
% Window = [-2 2];
% FigNum = 1;
% viewcell2b(cellid,'TriggerName','ResponseStart','SortEvent','ResponseEnd','eventtype','behav','ShowEvents',{[{'StimulusOffset'},{'ResponseEnd'}]},'ShowEventsColor',{[{'r'}]},'window',[-8 2],'PlotZeroLine','on','Partitions','#ChosenDirection');
%%
[psth spsth spsth_se tags spt stats] = ultimate_psth(cellid,'trial',TriggerEvent,[-2 4],'dt',0.01,'sigma',0.02,'parts','#ChosenDirection');
%%
xpsth1=-2:0.01:4;

figure
plot(xpsth1,spsth,'LineWidth',2)
set(gca,'FontSize',30,'box','off')
l1=legend('Left Choice','Right Choice');
set(l1,'FontSize',30,'box','off')
xlabel('Time from side port entry (s)','FontSize',25)
ylabel('Spike rate (sp/s)','FontSize',30)

%%
% Look at all the cells
for i =1:length(allcells)
cellid=allcells{i};

TE = loadcb(cellid,'TrialEvents');  % Behavioral data
SP = loadcb(cellid,'EVENTSPIKES'); % Spiking data

% Lets plot a psth


[TRIALS, TAGS, RTAGS, NUM_TAGS, PartNum] = partition_trials(TE,'#SideReward');

% Trigger on start of anticipation period
TriggerEvent = 'ResponseStart';

% SortingEvent='ResponseEnd';
% Window = [-2 2];
% FigNum = 1;
% viewcell2b(cellid,'TriggerName','ResponseStart','SortEvent','ResponseEnd','eventtype','behav','ShowEvents',{[{'StimulusOffset'},{'ResponseEnd'}]},'ShowEventsColor',{[{'r'}]},'window',[-8 2],'PlotZeroLine','on','Partitions','#ChosenDirection');

[psth spsth spsth_se tags spt stats] = ultimate_psth(cellid,'trial',TriggerEvent,[-2 4],'dt',0.01,'sigma',0.02,'parts','#ChosenDirection');

xpsth1=-2:0.01:4;

figure
plot(xpsth1,spsth,'LineWidth',2)
set(gca,'FontSize',30,'box','off')
l1=legend('Left Choice','Right Choice');
set(l1,'FontSize',30,'box','off')
xlabel('Time from side port entry (s)','FontSize',25)
ylabel('Spike rate (sp/s)','FontSize',30)
title(cellid,'FontSize',30)

end


%% Use the analysis to create a selection measure

% Side selectivity in the movement phase
EpochName='InitialResponse1';
nboot=200;


addanalysis(@calc_selectivity_sidePM,'property_names',{'Dside_Ant1';'Pside_Ant1'},'arglist',{EpochName;nboot});

%Outcome selectivity in olfactory trials in anticipation period
EpochName='InitialResponse1';
nboot=200;
addanalysis(@calc_selectivity_OutcOlfactoryPM,'property_names',{'Doutolf_Init1';'Poutolf_init1'},'arglist',{EpochName;nboot});

%Outcome selectivity in olfactory trials in anticipation period
EpochName='InitialResponse1';
nboot=200;
addanalysis(@calc_selectivity_OutcAuditoryPM,'property_names',{'Doutaud_Init1';'Poutaud_init1'},'arglist',{EpochName;nboot});


%% Look at Cell 1 

cellid=allcells{1};

TE = loadcb(cellid,'TrialEvents');  % Behavioral data
SP = loadcb(cellid,'EVENTSPIKES'); % Spiking data

% Lets plot a psth

% Trigger on start of anticipation period
TriggerEvent = 'ResponseStart';

% SortingEvent='ResponseEnd';
% Window = [-2 2];
% FigNum = 1;
% viewcell2b(cellid,'TriggerName','ResponseStart','SortEvent','ResponseEnd','eventtype','behav','ShowEvents',{[{'StimulusOffset'},{'ResponseEnd'}]},'ShowEventsColor',{[{'r'}]},'window',[-8 2],'PlotZeroLine','on','Partitions','#ChosenDirection');

[psth spsth spsth_se tags spt stats] = ultimate_psth(cellid,'trial',TriggerEvent,[-1 2],'dt',0.01,'sigma',0.02,'parts','#CorrectChoice');

xpsth1=-1:0.01:2;

figure
plot(xpsth1,spsth,'LineWidth',2)
set(gca,'FontSize',30,'box','off')
l1=legend('Error','Correct');
set(l1,'FontSize',30,'box','off')
xlabel('Time from side port entry (s)','FontSize',25)
ylabel('Spike rate (sp/s)','FontSize',30)
title(cellid,'FontSize',30)

%% Look at rates in this epoch

AntRates=SP.epoch_rates{2}; % Select teh rates from teh epoch

FullWaitTrials=TE.CompletedWTTrial==1 & TE.WaitingTime>2.8; % Select the trials


[FRFits,FRgof]=fit(AntRates(FullWaitTrials)',TE.WaitingTime(FullWaitTrials)','poly1');

xfr=linspace(min(AntRates(FullWaitTrials)),max(AntRates(FullWaitTrials)),1000);

yfr=FRFits.p1.*xfr+FRFits.p2;

cifr=predint(FRFits,xfr,0.95,'functional','on');

[r,p]=corrcoef(AntRates(FullWaitTrials),TE.WaitingTime(FullWaitTrials));
r(1,2)
p(1,2)

figure
scatter(AntRates(FullWaitTrials),TE.WaitingTime(FullWaitTrials),20,'k','LineWidth',1.4)
hold on
xlabel('Spike rate (sp/s)','FontSize',25)
ylabel('Waiting Time','FontSize',30)
title(cellid,'FontSize',30)
set(gca,'FontSize',30,'box','off')
plot(xfr,yfr,'b','LineWidth',2)
plot(xfr',cifr(:,1),'b')
plot(xfr',cifr(:,2),'b')



