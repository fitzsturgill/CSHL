% testing probe site screening for light effect....

Fs = 32000;
[fname, pname] = uiputfile('path', 'Choose CSC path...');
filepath = pname;

nlxcsc2mat2(filepath,'Channels','Events')
load(fullfile(filepath, 'Events.mat')); 

% event timestamps are in seconds, CSC timestamps are in microseconds
startRecording = find(strcmp(Events_EventStrings, 'Starting Recording'), 1, 'last');
startRecordingTime = Events_TimeStamps(startRecording);
stopRecording = find(strcmp(Events_EventStrings, 'Stopping Recording'), 1, 'last');
stopRecordingTime = Events_TimeStamps(startRecording);
% how many microseconds after recording start?



% h = waitbar(0, 'converting CSCs into TTs');   
spikeData = struct(...
    'spikes', [],...
    'amplitudes', []...
    );
spikeData = repmat(spikeData, 8, 1);

first = true;
% for trodeCounter = 1:8
nValidSamples = Nlx2MatCSC(fullfile(filepath, 'CSC1.ncs'),[0 0 0 1 0],0,1);%, [startRecording TrialStart_nlx(1)]);        
nValidSamples = sum(nValidSamples);
trode = 2;
ttch1 = 4 * (trode - 1) + 1;
chunkSize = 1e7;
nChunks = ceil(rem
    for sampleRange = 
    for ttcounter = 1:4
        cfilename = sprintf('CSC%d.ncs', ttch1 + ttcounter - 1);
%         [Samples, header] = Nlx2MatCSC(fullfile(filepath, cfilename),[0 0 0 0 1],1,1);%, [startRecording TrialStart_nlx(1)]);
        [nValidSamples, Samples, header] = Nlx2MatCSC(fullfile(filepath, cfilename),[0 0 0 1 1],1,1);%, [startRecording TrialStart_nlx(1)]);        
        display(['loaded ' cfilename]);
        sizeSamples = size(Samples);
        Samples = reshape(Samples, numel(Samples), 1);
        if first
            nSamples = length(Samples);
            first = false;
            ADBitVolts = sscanf(header{16}, '-ADBitVolts %f');            
            data = zeros(1, nSamples, 4);            
        end
        Samples = Samples * ADBitVolts * 1e6 * -1;
        
    %     % downsample to 20KHz
    %     Samples = resample(Samples, 32, 20);

        passBand = [300 8000]; % Hz
        [b, a] = butter(5, passBand/32000 * 2);
        Samples = filtfilt(b,a,Samples);
        data(1, :,ttcounter) = Samples;
    end
%     clear Samples
    display('detecting spikes');
%%    
    spikes = ss_default_params(Fs);
%%    
    spikes = ss_detect(data, spikes);
%%    
    spikeData(trode).spikes = ss_align(spikes);
%%    
    % compute spike amplitudes
    select = get_spike_indices(spikeData(trode).spikes, 'all');
    memberwaves = spikeData(trode).spikes.waveforms(select,:)';
    spikeData(trode).amplitudes = range(memberwaves);
    waitbar(trode/8);
% end
close(h);

%% plot amplitudes vs time
figname = 'spike_amplitudes_vs_time';
ensureFigure(figname, 1);
params.figmargin = [0.05 0.01 0.05 0.05];
fitz = colormap('jet');
nColors = size(fitz, 1);
colorIx = (1:8) * nColors / 8;
for trode = 1:8    
    axesmatrix(4,2,trode, params);    
    plot(spikeData(trode).spikes.spiketimes, spikeData(trode).amplitudes, 'Color', fitz(colorIx(trode), :)); 
    textBox(['tt' num2str(trode)], gca, [0.5 0.95], 12);
    set(gca, 'YLim', [0 300]);
    if ~rem(trode, 2)
        set(gca, 'YTick', []);
%         set(gca, 'YTickLabel', {''});
    end
    if trode < 7
        set(gca, 'XTick', []);
    end
    set(gca, 'XLim', [0 max(spikeData(trode).spikes.spiketimes)]);
    if trode == 1
        ylabel('spike amplitude (uV)');
    end
    if trode == 7
        xlabel('time (s)');
    end
end
    set(gcf, 'Position', [2162         335         753         602]);
saveas(gcf, fullfile(filepath, figname), 'jpeg');
saveas(gcf, fullfile(filepath, figname), 'fig');

%% throwaway, histograms of amplitudes and itis for 2 tetrodes, 

figname = 'spike_histograms';
ensureFigure(figname, 1);
subplot(1,2,1)
histogram(spikeData(6).amplitudes, 'FaceColor', fitz(colorIx(6), :)); hold on;
histogram(spikeData(1).amplitudes, 'FaceColor', fitz(colorIx(1), :)); 

set(gca, 'XLim', [0 300]);
ylabel('amplitude (uV)');
subplot(1,2,2)
histogram(diff(spikeData(1).spikes.spiketimes), 'FaceColor', fitz(colorIx(1), :)); hold on;
histogram(diff(spikeData(6).spikes.spiketimes), 'FaceColor', fitz(colorIx(6), :));
set(gca, 'XLim', [0 1]);
ylabel('ISI (s)');
saveas(gcf, fullfile(filepath, figname), 'jpeg');
saveas(gcf, fullfile(filepath, figname), 'fig');