function mcFindPeaks(fileCounter, automode, threshold)
% find peaks (i.e. for Multi-unit activity)
% automode == 0,     manual threshold, nargin must == 3 (threshold must be
% specified) 
% automode == 1,     autothreshold determined for each channel, 3xSD
% automode == 2,     autothreshold determined across all channels, 3xSD
    global state
    
    
    if nargin == 0 || isempty(fileCounter)
        fileCounter = 1:state.mcViewer.tsNumberOfFiles; % iterate through every loaded file
    end
    
    if nargin == 0
        automode = 1;
    end
    
    % determine automatic thresholds = 3x Std Dev
    % thresholds corresponds to [Matrix] of size (nChannels, nFilesToFindPeaks)
    if automode > 0 || nargin == 0
        invert = -1; %as default, look for negative goding peaks
        avgSD = mean(state.mcViewer.analysis.tsSD(:, fileCounter), 2);
        if automode == 2
            avgSD(:, 1) = mean(avgSD);
        end
        avgSD = repmat(avgSD, 1, length(fileCounter));
        state.mcViewer.analysis.tsAutoThresholds(:, fileCounter) =...
            state.mcViewer.analysis.tsMeans(:, fileCounter) + avgSD.*3;
        thresholds = state.mcViewer.analysis.tsAutoThresholds(:, fileCounter);
    else
        if size(threshold, 1) == 1
            thresholds = zeros(state.mcViewer.nChannels, length(fileCounter)) + threshold;
        else
            thresholds = repmat(threshold, 1, length(fileCounter));
        end
        
        if threshold > 0
            invert = 1;
        else
            invert = -1;
            thresholds = thresholds .* -1;
        end
    end
    

%     Do I even care to find the exact spike times at this point- i.e. the
%     average of rising and falling phases?  I can implement this later...
    lengthData = size(state.mcViewer.tsFilteredData{1,1}, 1);
    %find spike times and output times to tsSpikeTimes
    for i = 1:length(fileCounter)
        fileN = fileCounter(1, i);
        data = state.mcViewer.tsFilteredData{1, fileN} .* invert;
%         xData = repmat(state.mcViewer.tsXData{1, fileN}', 1,
%         state.mcViewer.nChannels);   
        for j = 1:state.mcViewer.nChannels
%               crossUpIndices = diff(state.mcViewer.tsFilteredData{1, fileN}(:, j) > thresholds(j, fileN)) == 1;
              crossUpIndices = diff(data(:, j) > thresholds(j, fileN)) == 1;
%         thisThresh = repmat(thresholds(:, i), 1, lengthData);
%             crossUpIndices = find(diff(data > thisThresh, 2) > 0);
%         crossDownIndices = find(diff(data > thisThresh, 2) < 0);
        
            crossUp =  state.mcViewer.tsXData{1, fileN}(crossUpIndices);
            state.mcViewer.analysis.tsSpikeTimes{j, fileN} = crossUp;
            state.mcViewer.analysis.tsSpikeAmps{j, fileN} = crossUp .* 0 + thresholds(j, fileN) .*invert;
            state.mcViewer.analysis.tsThresholds(j, fileN) = thresholds(j, fileN) .*invert;
        end
    end

%         crossDown = xData(crossDownIndices);
        
%         if any(min(crossDown, 1) < min(crossUp, 1));  % if there is one offending truncated spike, just deleted the first spike in all the traces- wrong!!!
%             crossUp = crossUp(2:end,:);
%             crossDown = crossUp(2:end-1, :);
%         end








%         SD = 0;
%         means = zeros(state.mcViewer.nChannels, length(fileCounter));  % divide
%         thresholds = zeros(state.mcViewer.nChannels, length(fileCounter));
%         %find average Std Dev across all files for each channel but means
%         %on a file by file basis
%         for i = 1:length(fileCounter)
%             fileN = fileCounter(1, i);
%             SD = SD + std(state.mcViewer.tsFilteredData{1, fileN})';
%             means(:,i) = mean(state.mcViewer.tsFilteredData{1, fileN})'; %now I've transposed
%             %does the previous line need to be transposed????
%         end
%         SD = SD./length(fileCounter); % avg SD across files for each channel
%         SD = repmat(SD, 1, length(fileCounter)); %extend SD dimensions along second dimension to emable matrix operation below
%         thresholds = means + SD.*3; % threshold = 3x STD DEV
        


        
        
        
        
    
    
    
        
    
    