nlxcsc2mat2(fullpth,'Channels',1:32);

%% CSC S/N determination...
animalID = 'CD4';
sessionID = '170912a';
fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];



%%
fullpth = 'Z:\tetrodeData\CD_3\2017-08-18_17-00-59';
rmsNoise = zeros(32, 1);
h = waitbar(0, 'RMS noise');   
for channel = 1:32
    cfilename = sprintf('CSC%d.ncs', channel);
    [Samples header] = Nlx2MatCSC(fullfile(fullpth, cfilename),[0 0 0 0 1],1,2,[1 1875]);
    Samples = reshape(Samples, numel(Samples), 1);s
    ADBitVolts = sscanf(header{16}, '-ADBitVolts %f');
    Samples = Samples * ADBitVolts * 1e6;
    %
    Fs = 32000;
%     window = [100 130];
%     i1 = bpX2pnt(window(1), Fs);
%     i2 = bpX2pnt(window(2), Fs);
%     xData = linspace(window(1), window(2), i2 - i1 + 1);
%     ensureFigure('cscData', 1); 
%     subplot(2,1,1);
%     plot(xData, Samples(i1:i2))
    filtData = bpFilterButterworth(Samples, 10000, 10, Fs);
%     subplot(2,1,2);
%     plot(xData, filtData);

    rmsNoise(channel) = rms(filtData);
    waitbar(channel/32);
end
close(h);
%%
ensureFigure('RMSNoise', 1);
plot(rmsNoise);




