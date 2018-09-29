%% Combine the scripts generating the files you need for further analysis
%% Add into cellbase, create waveform, etc
function AddSessionTagging_simple(Directory)

if nargin < 1
    Directory = uigetdir();
end
Directory = fullfile([Directory filesep]); % make sure it has a single trailing filesep

cd(Directory);

%% get start of recording time
nlxcsc2mat2(Directory,'Channels','Events')
load(fullfile(Directory, 'Events.mat'));
% 
% % event timestamps are in seconds, CSC timestamps are in microseconds
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

%% Creates and TT files
disp('Creating Waveform and TT files');
% Files=getDir(Directory,'file','clusters');
listing = dir;
Files = {};
for counter = 1:length(listing)
    if ~listing(counter).isdir && startsWith(listing(counter).name, 'clusters')
        Files{end + 1} = listing(counter).name;
    end
end
for f=1:length(Files)
    
    tetnum=sscanf(Files{f}, 'clusters%d.mat');
    
    ClusterFileName=strcat('clusters',num2str(tetnum),'.mat');
    CuratedFiringFile=strcat('firings',num2str(tetnum),'.curated.mda');
    
    load(ClusterFileName);
    MDAfile=readmda(CuratedFiringFile);
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
            
            
            save(TetrodeName,'tSpikes');
            
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
