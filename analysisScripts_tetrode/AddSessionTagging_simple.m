%% Combine the scripts generating the files you need for further analysis
%% Add into cellbase, create waveform, etc
function AddSessionTagging_simple(varargin)

%% optional parameters, first set defaults
defaults = {...
    'msdir', '';... % directory of mountainsort results (containing NT1,2,3.... subdirectories for each trode)
    'cbdir', '';... % cellbase session directory, i.e. \animal\session\  containing Events.nev file
    };
[s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings

if isempty(s.msdir)
    [~, msdir] = uiputfile('path', 'choose mountainsort directory, containing NT1,2,3... subdiretcories for each trode');
    if msdir == 0
        return
    else
        s.msdir = msdir;
    end
end

if isempty(s.cbdir)
    [~, cbdir] = uiputfile('path', 'choose cellbase session directory containing Events.nev file');
    if cbdir == 0
        return
    else
        s.cbdir = cbdir;
    end
end    


parts = strsplit(s.cbdir, filesep);
animal = parts{end - 1};
session = parts{end};

%% get start of recording time
nlxcsc2mat2(s.cbdir,'Channels','Events')
load(fullfile(s.cbdir, 'Events.mat'));

% event timestamps are in seconds, CSC timestamps are in microseconds
startRecording = find(strcmp(Events_EventStrings, 'Starting Recording'), 1, 'first');
startRecording = Events_TimeStamps(startRecording);



%% Extract the TTLS and create the Event and Event TTL Files
% Extract TTLs from neuralynx file
% events = getRawTTLs(fullfile(Directory,'Events.Nev'));

% EventTTL=events(:,2)';
% EventTimestamps=events(:,1)';
% EventTimestamps = EventTimestamps*10^-6;

% a=EventTTL==0;
% Array=find(a==1);
% Array=Array(1:end-1);
% foll=EventTTL(Array+1)==2;
% 
% EventTTL(Array(foll)+1)=999;
% 
% EventTTL(EventTTL==2)=998;
% EventTTL(EventTTL==999)=2;
% 
% save('EventTTL_task.mat','EventTTL');
% save('Events.mat','EventTimestamps','EventTTL');

%% Creates TT files in cellbase directory
disp('Creating Waveform and TT files');


cd(s.msdir);
trodeFolders = dir(['**' filesep 'NT*']); % dir('**\NT*')-  find NT1,2,3.... folders for each trode


for counter=1:length(trodeFolders)
    if ~trodeFolders(counter).isdir
        continue % this should never happen
    end
    
    tetnum=sscanf(trodeFolders(counter).name, 'NT%d');
    
%     ClusterFileName=strcat('clusters',num2str(tetnum),'.mat');
%     CuratedFiringFile=strcat('firings',num2str(tetnum),'.curated.mda');
    
    load(fullfile(trodeFolders(counter).folder, trodeFolders(counter).name, filesep, 'clusters.mat'));
    try
        MDAfile=readmda(fullfile(trodeFolders(counter).folder, trodeFolders(counter).name, filesep, 'firings.curated.mda'));
    catch
        continue;
    end
    
    if ~isempty(MDAfile)
        cells=unique(MDAfile(3,:));
        numcells=length(cells);
        
        
        for nc=1:numcells
            Index=(MDAfile(3,:)==cells(nc));
            Times=ismember(clusters.firings(2,:),MDAfile(2,Index));
%             TS=clusters.firings_seconds(Times);
            %             TS=TS*10000; % Change the timing depending on the time basis
            tSpikes = clusters.firings(2,Times);
            tSpikes = tSpikes / 32000; % convert samples into seconds given 32000Hz sample rate
            tSpikes = tSpikes + startRecording;
            TetrodeName = sprintf('%s%d_%02d.mat','TT',tetnum,nc);
            
            
            save(fullfile(s.cbdir, TetrodeName),'tSpikes');
            
        end
    end
end

disp('TT files created');

% %% Create stimulus events
% MakeStimEvents(Directory,'PulseNttl',526); %
% 
% [part,sess]=fileparts(Directory);
% [~,anim]=fileparts(part);
% SubDirectory=fullfile(anim,sess);
% CellList=findallcells(SubDirectory);
% 
% % Prealign spikes to stimulus events
% problem_stim_cellid = [];
% for iC = 1:length(CellList)
%     cellid = CellList(iC);
%     try
%         prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_laserstim,'filetype','stim','ifsave',1,'ifappend',0,'writing_behavior','overwrite')
%     catch
%         disp('Error in prealignSpikes.');
%         problem_stim_cellid = [problem_stim_cellid cellid];
%     end
% end
% 
% % View light-triggered raster and PSTH
% TrigEvent = 'PulseOn';
% SEvent = 'PulseOn';
% win = [-0.05 0.05];
% parts = 'all';
% %     parts = '#BurstNPulse';
% dt = 0.001;
% sigma = 0.001;
% PSTHstd = 'on';
% ShEvent = {'PulseOn'};
% ShEvColors = hsv(length(ShEvent{1}));
% ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
% for iCell = 1:length(CellList)
%     cellid = CellList{iCell};
%     H = figure;
%     viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
%         'FigureNum',1,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
%         'EventMarkerWidth',0,'PlotZeroLine','off');
% end
% 
% end
% 
