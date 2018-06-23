



[fname, pname] = uiputfile('path', 'Choose CSC path...');
filepath = pname;
nlxcsc2mat2(filepath,'Channels','Events')
load(fullfile(filepath, 'Events.mat'));

%%
% laser pulse times
PortID = eventPortFromEventStrings(Events_EventStrings);
Pulses = PortID == 0 & Events_Nttls == 128;
PulseTimes = Events_TimeStamps(Pulses);
% get a window of 10 pulses
timeStampRange = [PulseTimes(1) PulseTimes(100)] * 1e6;

%%
channels = 1:4:32;
nChannels = length(channels);
for counter = 1:length(channels)
    channel = channels(counter);
    cfilename = sprintf('CSC%d.ncs', channel);
    disp(cfilename);
%     TimeStamps, ChannelNumbers, SampleFrequencies, NumberValidSamples, Samples and NlxHeader 
%     [TimeStamps, Samples, header] = Nlx2MatCSC(fullfile(filepath, cfilename),[1 0 0 0 1],1,1);% testing all the data, 1 channel
    [TimeStamps, Samples, header] = Nlx2MatCSC(fullfile(filepath, cfilename),[1 0 0 0 1],1,4, timeStampRange);%, [startRecording TrialStart_nlx(1)]);
    Samples = reshape(Samples, numel(Samples), 1);
    ADBitVolts = sscanf(header{16}, '-ADBitVolts %f');
    Samples = Samples * ADBitVolts * 1e6;    
    if channel == 1
        data = zeros(numel(Samples), nChannels);
    end
    data(:,counter) = Samples;
end

%%
nSamples = size(data,1);
xdata = linspace(0, nSamples / 32000, nSamples);
ensureFigure('artifact', 1);
for counter = 1:8    
    subplot(3,3,counter);
    plot(xdata, data(:,counter));
    legend(['channel ' num2str(channels(counter))]);
%         set(gca, 'XLim', [0 0.1]);
end
set(gcf, 'Position', [518 107 1136 789]);
saveas(gcf, fullfile(filepath, 'artifact2'), 'jpeg');
saveas(gcf, fullfile(filepath, 'artifact2'), 'fig');

%%
ensureFigure('artifact_zoom', 1);
for counter = 1:8    
    subplot(3,3,counter);
    plot(xdata, data(:,counter));
    legend(['channel ' num2str(channels(counter))]);
    set(gca, 'XLim', [0.01 0.015]);
end
set(gcf, 'Position', [518 107 1136 789]);
saveas(gcf, fullfile(filepath, 'artifact_zoom'), 'jpeg');
saveas(gcf, fullfile(filepath, 'artifact_zoom'), 'fig');


