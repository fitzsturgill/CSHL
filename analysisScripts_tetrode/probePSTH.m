function spikeTimes = probePSTH(varargin)

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
        'avg_window', [-0.1 0.1];...  % 100msec snippet after laser pulse start to generate CSC avg per lead.
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
PulseTimes = Events_TimeStamps(Pulses);

h = waitbar(0, 'Detecting Spikes');   
nValidSamples = Nlx2MatCSC(fullfile(s.filepath, 'CSC1.ncs'),[0 0 0 1 0],0,1);% get # total sample pieces (output is nSamplePieces x 512)
% sampleCheck.nlx = sum(nValidSamples);
nValidSamples = numel(nValidSamples);
chunkSize = 1e12; % 1e5
nChunks = ceil(nValidSamples / chunkSize);
sampleRange = ((0:nChunks) * chunkSize);

%%
passBand = [300 8000]; % Hz
[z, p, k] = butter(5, passBand/(s.Fs/2));
[sos, g] = zp2sos(z,p,k);
spikeTimes = cell(s.nChannels, 1);
first = 1;

for channel = 1:s.nChannels
    spikes = ss_default_params(s.Fs, {'display', 'trial_spacing'}, 0); % don't pad time between chunks/trials, see ss_detect
    for counter = 1:nChunks
        theseSamples = [(sampleRange(counter) + 1) min(sampleRange(counter + 1), nValidSamples)];          
    %     sampleCheck.range(counter, :) = theseSamples;
            cfilename = sprintf('CSC%d.ncs', channel);
            [Timestamps, Samples, header] = Nlx2MatCSC(fullfile(s.filepath, cfilename),[1 0 0 0 1],1,2, theseSamples); %, [startRecording TrialStart_nlx(1)]);        
            %         find the time offset for the current chunk
%             offset = TimeStamps(1);                        
            %             display(['loaded ' cfilename]);
            % get length of data chunk (last chunk is shorter), bitVolts
            % conversion factor, initialize data array for spike detection
            Samples = reshape(Samples, numel(Samples), 1);        
            if first                
                ADBitVolts = sscanf(header{16}, '-ADBitVolts %f');                            
            end
            Samples = Samples * ADBitVolts * 1e6 * inv;
            Samples = filtfilt(sos,g,Samples);            
            spikes = ss_detect({Samples}, spikes);            
    end
%         select = get_spike_indices(spikes, show );     
%     memberwaves = spikes.waveforms(select,:);
%     spiketimes  = sort( spikes.unwrapped_times( select ) );
%     % get amplitudes over time
%     amp = range(    memberwaves' );
    spikeTimes{channel,1} = double(spikes.spiketimes(:)) + double(startRecording);           

    sprintf('processed channel %d\n', channel)
    waitbar(channel/s.nChannels);
    break;
end

close(h);




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

