function stimData = probePSTH2(varargin)

    % generate probe laser PSTH from each lead, also compute transformation matrix for ZCA
    % whitening
    % for laser PSTH- subtract mean of signal across all the other leads
    % prior to spike detection
    % also make average CSC for all the laser pulses....
    
    defaults = {...
        'direction', 'up';...
        'nChannels', 32;...
        'Fs', 32000;...
        'filepath', '';...
        'avg_window', [-0.1 0.1];...  % if this extends left and/or right beyond adjacent pulses in a burst, it will merge all the pulses together into a "trial" segement for spike sorting
        };
    
    [s, ~] = parse_args(defaults, varargin{:});
    
    if isempty(s.filepath)
        s.filepath = uigetdir();
        if isempty(s.filepath)
            return
        end
    end
    s.filepath = fullfile([s.filepath filesep]);    
    switch s.direction
        case 'up'
            inv = -1;
        case 'down'
            inv = 1;
    end
    
% testing probe site screening for light effect....
nlxcsc2mat2(s.filepath,'Channels','Events');
load(fullfile(s.filepath, 'Events.mat'));

% % event timestamps are in seconds, CSC timestamps are in microseconds
startRecording = find(strcmp(Events_EventStrings, 'Starting Recording'), 1, 'first');
startRecording = Events_TimeStamps(startRecording);
TrialStart_nlx = getBehaviorStartTimes(Events_Nttls, Events_EventStrings, Events_TimeStamps);

PortID = eventPortFromEventStrings(Events_EventStrings);
Pulses = PortID == 0 & Events_Nttls == 128;
if sum(Pulses) == 0
    disp(sprintf('no pulses found for %s', fullfile(s.filepath, 'Events.mat')));
    stimData = [];
    return
end
PulseTimes = Events_TimeStamps(Pulses);

% get scaling factor
header = Nlx2MatCSC(fullfile(s.filepath, 'CSC1.ncs'),[0 0 0 0 0],1,1);
ADBitVolts = sscanf(header{16}, '-ADBitVolts %f');


% build up spike detection segments (merging overlapping windows) that comprise portions of
% recording containing laser pulses for spike detection

segments = [];
oldWindow = PulseTimes(1) + s.avg_window;
nPulses = length(PulseTimes);
for counter = 2:length(PulseTimes)
    thisWindow = PulseTimes(counter) + s.avg_window;    
    
    if thisWindow(1) > oldWindow(2)
        segments(end + 1, :) = oldWindow;
        oldWindow = thisWindow;
    else
        oldWindow(2) = thisWindow(2);
    end
    
    if counter == nPulses
        segments(end + 1, :) = oldWindow; % add the last one no matter what
    end
end

h = waitbar(0, 'Detecting Spikes');   
ValidSamples = Nlx2MatCSC(fullfile(s.filepath, 'CSC1.ncs'),[0 0 0 1 0],0,1);% get # total sample pieces (output is nSamplePieces x 512)
% sampleCheck.nlx = sum(nValidSamples);
nValidSamples = sum(ValidSamples);
% chunkSize = 1e12; % 1e5
% nChunks = ceil(nValidSamples / chunkSize);
% sampleRange = ((0:nChunks) * chunkSize);

%%
passBand = [300 8000]; % Hz
[z, p, k] = butter(5, passBand/(s.Fs/2));
[sos, g] = zp2sos(z,p,k);

% first = 1;

% data to output, initialize structure:
% spike numbers, spike amplitudes, CSC snippet per light pulse for
% subsequent averaging,     also, option to gather random subset for ZCA
% whitening


stimData = struct(...
    'spikeTimes', [],...
    'spikeAmplitudes', [],...
    'snippet', []... % try downsampling 10x and storing for every channel
    );
stimData = repmat(stimData, 32, 1);

nSegments = size(segments, 1);
subTimeStamps = 0:1/s.Fs:(1/s.Fs * 511);
subTimeStamps = subTimeStamps(:);

% for P-2 probe using neurolynx adapter
channelOrder = [16 5 12 4 10 2 9 1 7 3 8 6 14 13 15 11 17 28 21 29 23 31 24 32 26 30 25 27 19 20 18 22];

% for E-2 probe:
channelOrder = [16 12 10 9 7 5 4 2 15 11 14 13 8 6 3 1 17 21 23 24 26 28 29 31 18 22 19 20 25 27 30 32];
for channel = 1:s.nChannels
    tic
    cfilename = sprintf('CSC%d.ncs', channelOrder(channel)); 
    cscData = cell(nSegments, 2); % first dimension is voltage (Samples), second dimension is time (Timesteps)
    downData = []; % downsampled for making averages later
    trialStarts = zeros(nSegments, 1);
    spikes = ss_default_params(s.Fs, {'display', 'trial_spacing'}, 0); % don't pad time between chunks/trials, see ss_detect
    for counter = 1:nSegments
        [Timestamps, Samples, header] = Nlx2MatCSC(fullfile(s.filepath, cfilename),[1 0 0 0 1],1,4, segments(counter,1:2) * 1e6); % convert event timestamps to microseconds 
        downData = expandVertCat(downData, decimate(Samples(:) * ADBitVolts * 1e6 * inv, 10)');
        Samples = filtfilt(sos,g,Samples(:) * ADBitVolts * 1e6 * inv);
        cscData{counter, 1} = Samples;
        Timestamps_full = repmat(Timestamps, 512, 1) + repmat(subTimeStamps, 1, numel(Timestamps));
        cscData{counter, 2} = Timestamps_full(:);
        trialStarts(counter) = Timestamps(1);
    end
    spikes = ss_detect(cscData(:,1), spikes);
    toc
    stimData(channel).spikeTimes = spikes.spiketimes;
    stimData(channel).spikeAmplitudes = range(spikes.waveforms');
    stimData(channel).snippet = downData;
    sprintf('processed channel %d\n', channelOrder(channel))
    waitbar(channel/s.nChannels);
    
    %     memberwaves = spikes.waveforms(select,:);
%     spiketimes  = sort( spikes.unwrapped_times( select ) );
%     % get amplitudes over time
%     amp = range(    memberwaves' );
% if channel == 10
%     break
% end
end
% stuff.segments = segments;
% stuff.spikes = spikes;
% stuff.cscData = cscData;
% stuff.trialStarts = trialStarts;
% stuff.PulseTimes = PulseTimes;
close(h);


%     spikeTimes{channel,1} = double(spikes.spiketimes(:)) + double(startRecording);           

    save(fullfile(s.filepath, 'stimData.mat'), 'stimData');
    disp(['*** Saved: ' fullfile(s.filepath, 'stimData.mat')]);


return;
% %% scrapbook for converting into t files for cellbase
%% for each channel, calculate laser PSTH
nChannels = 32;
PortID = eventPortFromEventStrings(Events_EventStrings);
Pulses = PortID == 0 & Events_Nttls == 128;
PulseTimes = Events_TimeStamps(Pulses);

% 10Hz stimulation, 10ms bins
edges  = -0.05:.01:0.05;
histCountsByChannel = zeros(length(edges) - 1, nChannels);

for counter = 1:length(PulseTimes)
    pulseTime = PulseTimes(counter);
    for channel = 1:nChannels
        spikesZeroed = spikeTimes{channel} - pulseTime;
        histCountsByChannel(:,channel) = histCountsByChannel(:,channel) + histcounts(spikesZeroed, edges)';        
    end
end

histCountsByChannel_norm = (histCountsByChannel - mean(histCountsByChannel(1:3, :), 1)) ./ mean(histCountsByChannel(1:3, :), 1);

