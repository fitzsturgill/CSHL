[fname, pname] = uiputfile('path', 'Choose CSC path...');
filepath = pname;



nChannels = 32;


data = zeros(32000 * 10, 32);
passBand = [300 8000]; % Hz
[b, a] = butter(5, passBand/32000 * 2);
for channel = 1:nChannels
    disp(['loading channel ' num2str(channel)]);
    cfilename = sprintf('CSC%d.ncs', channel);
%     TimeStamps, ChannelNumbers, SampleFrequencies, NumberValidSamples, Samples and NlxHeader 
%     [TimeStamps, Samples, header] = Nlx2MatCSC(fullfile(filepath, cfilename),[1 0 0 0 1],1,2, [1 1875 * 2]);% 20s of data
%     [a, b, c, d, e] = Nlx2MatCSC(fullfile(filepath, cfilename),[1 1 1 1 1],1,2, [1 1875 * 2]);% 20s of data
%     [TimeStamps, Samples, header] = Nlx2MatCSC(fullfile(filepath, cfilename),[1 0 0 0 1],1,1);% testing all the data, 1 channel
[NumValSamples, header, Samples] = Nlx2MatCSC(fullfile(filepath, cfilename),[0 0 0 1 1],1,1);% testing all the data, 1 channel
    break
    Samples = reshape(Samples, numel(Samples), 1);
    ADBitVolts = sscanf(header{16}, '-ADBitVolts %f');
    Samples = Samples * ADBitVolts * 1e6;

    filtData = filtfilt(b,a,Samples);
    if channel == 1
        data = zeros(numel(filtData), 32);
    end
    data(:,channel) = filtData;    
end   

%% 
% [Xwh, mu, Winv, W] = whiten(data);
% 
% %%
% ensureFigure('testWhiten', 1);
% channel = 14;
% drange = [1 32000] + 10 * 32000;
% 
% % for counter = 1:nRows
% 
% %     subplot(2, 1, 1);
%     plot(zscore(data(drange(1):drange(2), channel)));
% %     subplot(2, 1, 2); 
% hold on;
%     plot(zscore(Xwh(drange(1):drange(2), channel)));
% % end

    
