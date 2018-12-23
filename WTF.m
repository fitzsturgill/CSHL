basepath = 'C:\Users\Fitz\Documents\Data_Local\iChR2_test\';

folders = dir(basepath);
folders = folders(3:end);


for counter = 2:length(folders)
%     try
        probePSTH2('filepath', fullfile(basepath, folders(counter).name, filesep), 'avg_window', [-3 4]);
%     catch ME
        disp(fullfile(basepath, folders(counter).name, filesep));
%     end
end
    
%%
return;
%%
ensureFigure('test', 1); 

offset = 0;
for counter = 1:32
    subplot(4,8,counter);
    scatter(stimData(counter + offset).spikeTimes, stimData(counter + offset).spikeAmplitudes, '.');
    textBox(['ch' num2str(counter + offset)]);
end


%%

ensureFigure('test2', 1); 

offset = 0;
for counter = 1:32
    subplot(16,2,counter);
    plot(linspace(-1, 2, size(stimData(counter).snippet, 2)), mean(stimData(counter).snippet));
    textBox(['ch' num2str(counter + offset)]);
%     set(gca, 'XLim', [3.2 4.2]); set(gca, 'YLim', [-200 200]);
end
