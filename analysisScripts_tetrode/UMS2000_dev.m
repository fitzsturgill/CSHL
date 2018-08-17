function [spikes sampleCheck] = UMS2000_dev(direction)
    global sampleCheck
    if nargin < 1
        direction = 'up';
    end
    
    switch direction
        case 'up'
            inv = -1;
        case 'down'
            inv = 1;
    end
% testing probe site screening for light effect....

Fs = 32000;
% [fname, pname] = uiputfile('path', 'Choose CSC path...');
% filepath = pname;
filepath = 'F:\Cellbase\CD17\180710a\';

nlxcsc2mat2(filepath,'Channels','Events')
load(fullfile(filepath, 'Events.mat'));
% 
% % event timestamps are in seconds, CSC timestamps are in microseconds
startRecording = find(strcmp(Events_EventStrings, 'Starting Recording'), 1, 'first');
startRecording = Events_TimeStamps(startRecording);
% stopRecording = find(strcmp(Events_EventStrings, 'Stopping Recording'), 1, 'last');
% stopRecordingTime = Events_TimeStamps(startRecording);
% how many microseconds after recording start?



h = waitbar(0, 'Detecting Spikes');   
% spikeData = struct(...
%     'spikes', [],...
%     'amplitudes', []...
%     );
% spikeData = repmat(spikeData, 8, 1);


nValidSamples = Nlx2MatCSC(fullfile(filepath, 'CSC1.ncs'),[0 0 0 1 0],0,1);% get # total sample pieces (output is nSamplePieces x 512)
sampleCheck.nlx = sum(nValidSamples);
nValidSamples = numel(nValidSamples);
chunkSize = 1e5;
nChunks = ceil(nValidSamples / chunkSize);
sampleRange = ((0:nChunks) * chunkSize);
% for trodeCounter = 1:8

trode = 2;
ttch1 = (4 * (trode - 1)) + 1;

sampleCheck.perChunk = zeros(nChunks, 1);
sampleCheck.range = zeros(nChunks, 2);
spikes = ss_default_params(Fs, {'display', 'trial_spacing'}, 0); % don't pad time between chunks/trials, see ss_detect
%% debugging
% segSize = 100;
% sampleCheck.beforeSegments = zeros(nChunks, segSize);
% sampleCheck.afterSegments = sampleCheck.beforeSegments;
%%

for counter = 1:nChunks
    theseSamples = [(sampleRange(counter) + 1) min(sampleRange(counter + 1), nValidSamples)]
    sampleCheck.range(counter, :) = theseSamples;
    for ttcounter = 1:4
        cfilename = sprintf('CSC%d.ncs', ttch1 + ttcounter - 1);
        [Timestamps, Samples, header] = Nlx2MatCSC(fullfile(filepath, cfilename),[1 0 0 0 1],1,2, theseSamples); %, [startRecording TrialStart_nlx(1)]);        
        display(['loaded ' cfilename]);
        % get length of data chunk (last chunk is shorter), bitVolts
        % conversion factor, initialize data array for spike detection
        Samples = reshape(Samples, numel(Samples), 1);        
        if ttcounter == 1
            nSamples = length(Samples);
            ADBitVolts = sscanf(header{16}, '-ADBitVolts %f');            
            data = zeros(1, nSamples, 4);
        end
        sampleCheck.perChunk(counter) = numel(Samples);
        %% debugging
%         if ttcounter == 1
%             disp(cfilename);
%             sampleCheck.beforeSegments(counter, :) = Samples(1:segSize);
%             sampleCheck.afterSegments(counter, :) = Samples(end - segSize + 1:end);
%         end
        %%
        
        Samples = Samples * ADBitVolts * 1e6 * inv;

        
    %     % downsample to 20KHz
    %     Samples = resample(Samples, 32, 20);

        passBand = [300 8000]; % Hz
        [b, a] = butter(5, passBand/Fs * 2);
        Samples = filtfilt(b,a,Samples);
        data(1, :,ttcounter) = Samples;
        waitbar(((counter - 1) * 4 + ttcounter) / (nChunks * 4));
    end   
    spikes = ss_detect(data, spikes);
    disp('spikes detected');
end
% cfilename = sprintf('CSC%d.ncs', 5);
% [Samples, header] = Nlx2MatCSC(fullfile(filepath, cfilename),[0 0 0 0 1],1,1);
% Samples = reshape(Samples, numel(Samples), 1);
% sampleCheck.samples = Samples;

sampleCheck.chunksSummed = sum(sampleCheck.perChunk);
spikes.startRecording = startRecording;
close(h);
spikes = ss_align(spikes);
spikes = ss_kmeans(spikes);
spikes = ss_energy(spikes);
spikes = ss_aggregate(spikes);


%% scrapbook for converting into t files for cellbase
return
show = spikes.assigns == 155;
tSpikes = spikes.unwrapped_times(show);
tSpikes = double(tSpikes');
tSpikes = tSpikes + spikes.startRecording;

