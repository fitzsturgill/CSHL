% test_script for 
% loadcb; cellid=CELLIDLIST{1}; tagging(cellid); taggedpropTO(cellid); taggedsummaryTO(cellid)

% figure; viewcell2b(thiscell,'TriggerName','Us_start','SortEvent','trialNumber','eventtype','behav','ShowEvents',{'Us_start'},...
%    'Partitions','#trialType: {1 4 7}','window',[-7 4], 'dt', 0.01, 'sigma', 0.02, 'PSTHstd', 'on', 'isadaptive', true);

saveOn = 1;
%%
rat = 'CD21';
session = '181221a';
loadcb; % % loads ANALYSES, CELLIDLIST, PREFERENCES, and TheMatrix
pref = getcbpref; % torben's new way
fullpath = [pref.datapath '\' rat '\' session '\'];
%%
addnewcells('dir',[rat filesep session]);
nlxcsc2mat2(fullpath,'Channels','Events');
TE = makeTE_CuedOutcome_Odor_Complete_Nlx(fullpath);
MakeStimEvents_Bpod(fullpath,'PulseNttl',128, 'PulsePort', 0); % FS


%% add number of outcome Licks occuring between 0.1 and 1s following outcome
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 1], TE.Us);
TE.usLicks_rate = TE.usLicks.rate;
TE.csLicks = countEventFromTE(TE, 'Port1In', [0 3], TE.Cue);
TE.csLicks_rate = TE.csLicks.rate;
TE.rewardReceived = TE.usLicks.count >= 2; % get rid of trials where mouse fails to lick for reward
save(fullfile(pref.datapath, rat, session, pref.TrialEvents_fname), 'TE'); 

%%


thesecells = findcell('rat', rat, 'session', session);

%%
for counter = 1:length(thesecells)
    cellid = thesecells{counter};
    % [ratname,session,tetrode,unit] = cellid2tags(cellid);
    TE = loadcb(cellid, 'TrialEvents');

    % eventIx = find(strcmp(SP.events(:,1), 'Us_start'));
    % stimes = SP.event_stimes{eventIx};

    %
    prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_CuedOutcome,'filetype','event','ifsave',1,'ifappend',0, 'writing_behavior', 'overwrite');
    prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_laserstim_Bpod,'filetype','stim','ifsave',1,'ifappend',0);
    SP = loadcb(cellid,'EVENTSPIKES');
    % make rasters and PSTHs
    % time = window(1) - winMargin:dt:window(2) + winMargin;
    % binraster = stimes2binraster(stimes,time,dt);

    %
    cellidStripped = regexprep(cellid,'\.','_');
    window = [-7 4];
    dt = 0.01;
    sigma = 0.02;
    winMargin = sigma * 3;
    isadaptive = false;

    fh = [];
    figName = [cellidStripped '_rewardSpikes'];    
    ensureFigure(figName, 1);  viewcell2b(cellid,'TriggerName','Us_start','SortEvent','csLicks_rate',...
        'eventtype','behav','ShowEvents',{'Us_start', 'Cue_start'}, 'Partitions','#trialType: {1 4 7} & rewardReceived',...
        'window',[-2 2], 'dt', dt, 'sigma', sigma, 'PSTHstd', 'on', 'isadaptive', isadaptive);
    legend('High value', 'Low value', 'Uncued'); title(figName, 'Interpreter', 'none');
    fh(end+1) = gcf; formatFigureCellbase;

    figName = [cellidStripped '_punishSpikes'];
    ensureFigure(figName, 1);  viewcell2b(cellid,'TriggerName','Us_start','SortEvent','trialNumber',...
        'eventtype','behav','ShowEvents',{'Us_start', 'Cue_start'}, 'Partitions','#trialType: {2 5 8}',...
        'window',[-2 2], 'dt', dt, 'sigma', sigma, 'PSTHstd', 'on', 'isadaptive', isadaptive);
    legend('High value', 'Low value', 'Uncued'); title(figName, 'Interpreter', 'none');
    fh(end+1) = gcf; formatFigureCellbase;

    figName = [cellidStripped '_rewardLicks'];
    ensureFigure(figName, 1); viewlick({rat, session}, 'TriggerName', 'Us_start', 'SortEvent', 'csLicks_rate',...
        'eventtype', 'behav', 'ShowEvents', {'Us_start', 'Cue_start'}, 'Partitions', '#trialType: {1 4 7} & rewardReceived',...
        'window',[-2 2], 'dt', dt, 'sigma', sigma, 'PSTHstd', 'on', 'isadaptive', isadaptive, 'LickInField', 'Port1In');
    legend('High value', 'Low value', 'Uncued'); title(figName, 'Interpreter', 'none');
    fh(end+1) = gcf; formatFigureCellbase;

    figName = [cellidStripped '_cueSpikes'];
    ensureFigure(figName, 1);  viewcell2b(cellid,'TriggerName','Cue_start','SortEvent','csLicks_rate',...
        'eventtype','behav','ShowEvents',{'Cue_start'}, 'Partitions','#cueCondition',...
        'window',[-2 5], 'dt', dt, 'sigma', sigma, 'PSTHstd', 'on', 'isadaptive', isadaptive);
    title(figName, 'Interpreter', 'none');
    fh(end+1) = gcf; formatFigureCellbase;

    figName = [cellidStripped '_cueLicks'];
    ensureFigure(figName, 1); viewlick({rat, session}, 'TriggerName', 'Cue_start', 'SortEvent', 'csLicks_rate',...
        'eventtype', 'behav', 'ShowEvents', {'Cue_start'}, 'Partitions', '#cueCondition',...
        'window',[-2 5], 'dt', dt, 'sigma', sigma, 'PSTHstd', 'on', 'isadaptive', isadaptive, 'LickInField', 'Port1In');
    title(figName, 'Interpreter', 'none');
    fh(end+1) = gcf; formatFigureCellbase;

    figName = [cellidStripped '_outcomeSpikes'];    
    ensureFigure(figName, 1);  viewcell2b(cellid,'TriggerName','Us_start', 'SortEvent', 'trialNumber',...
        'eventtype','behav','ShowEvents',{'Us_start'}, 'Partitions','#trialOutcome',...
        'window',[-2 2], 'dt', dt, 'sigma', sigma, 'PSTHstd', 'on', 'isadaptive', isadaptive);
    title(figName, 'Interpreter', 'none');
    fh(end+1) = gcf; formatFigureCellbase;

    % tagging
    figName = [cellidStripped '_taggingSpikes'];   
    ensureFigure(figName, 1);  
    viewcell2b(cellid,'TriggerName','BurstOn','SortEvent','BurstOn','ShowEvents',{'PulseOn'},...
        'eventtype','stim','window',[-2 2],'dt',dt,'sigma',sigma,'PSTHstd','on',...
        'EventMarkerWidth',0,'PlotZeroLine','off', 'isadaptive', isadaptive);
    title(figName, 'Interpreter', 'none');
    fh(end+1) = gcf; formatFigureCellbase;


    pdfname = fullfile(fullpath,['GradedValueSummary_' regexprep(cellid,'\.','_') '.pdf']);
    if saveOn
        for counter = 1:length(fh)
            if counter == 1
                export_fig(fh(counter),pdfname);  % write to pdf
            else
                export_fig(fh(counter),'-append',pdfname);  % write to pdf
            end

        end
    end
end
%%



% MakeStimEvents_Bpod(fullpath,'PulseNttl',128, 'PulsePort', 0); % FS
% prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_CuedOutcome,'filetype','event','ifsave',1,'ifappend',0, 'writing_behavior', 'overwrite')


