function mcFilterData(files,force,channels)


    global state
    if nargin == 0
        files = state.mcViewer.fileCounter;
    end
    
    % if data is already filtered don't filter mcChannels unless force = 1
    if nargin < 2
        force = 0;
    end
    
    % if channels aren't specified
    if nargin < 3
        channels = 1:state.mcViewer.totalChannels;
    end
    
    % no data is loaded to filter
    if isempty(state.mcViewer.tsData)
        return
    end

    if length(files) > 1
        h = waitbar(0, 'Filtering Data');
    end
    for i = files
        %% Filter Data- first mcChannel Data

        % if any mcChannel is specified filter them all and continue to
        % auxiliary channels.  mcChannel filtering is controlled separately
        % from channelControl GUI elements.  ChannelConrtrol GUI elements only
        % determine display filtering for these channels.
        
        % only filter fileNumber if data is not already filtered or force is
        % specified
        if force || size(state.mcViewer.tsFilteredData, 2) < i || isempty(state.mcViewer.tsFilteredData{1, i})        
            if any(channels < state.mcViewer.nChannels)

                lowCutoff = 0;
                highCutoff = 0;
                if state.mcViewer.lowPass
                    lowCutoff = state.mcViewer.lowCutoff;
                end
                if state.mcViewer.highPass
                    highCutoff = state.mcViewer.highCutoff;
                end

                filtData = mcFilter(state.mcViewer.tsData{1, i}(:,1:state.mcViewer.nChannels), lowCutoff, highCutoff, state.mcViewer.sampleRate);
                state.mcViewer.analysis.tsSD(:,i) = std(filtData, 0, 1)';
                state.mcViewer.analysis.tsMeans(:, i) = mean(filtData, 1)';
                if state.mcViewer.blankTrain
                    pulseTrain = repmat(mcMakePulseTrain', 1, state.mcViewer.nChannels);
                    if isempty(pulseTrain)
                        state.mcViewer.tsFilteredData{1, i}(:, 1:state.mcViewer.nChannels) = filtData;
                    else
                        state.mcViewer.tsFilteredData{1, i}(:, 1:state.mcViewer.nChannels) = filtData .* pulseTrain; %FS MOD Kludge- need to add GUI controls for this...
                    end
                else
                    state.mcViewer.tsFilteredData{1, i}(:, 1:state.mcViewer.nChannels) = filtData;
                end
            end
            %% Filter Data- second, auxiliary channels, controlled by
            % ChannelControl GUI elements
            auxChannels = channels(channels > state.mcViewer.nChannels);
            for j = auxChannels
                if state.mcViewer.channel(j).ShowFilter
                    % filter the channel
                    state.mcViewer.tsFilteredData{1, i}(:, j) = ...
                        mcFilter(state.mcViewer.tsData{1, i}(:, j), ...
                        state.mcViewer.channel(j).LowPass, ...
                        state.mcViewer.channel(j).HighPass ...
                        );          
                else
                    state.mcViewer.tsFilteredData{1, i}(:, j) = ...
                                    state.mcViewer.tsData{1, i}(:, j);
                end
            end
        end
        if length(files) > 1
            mcFlipTimeSeries(i);
            waitbar(i/length(files));
        end
    end
    if length(files) > 1
        close(h);
    end
% function mcFilterData(fileCounter,force,channels)
% 
% 
%     global state
%     if nargin == 0
%         fileCounter = state.mcViewer.fileCounter;
%     end
%     
%     % if data is altready filtered return unless force = 1
%     if nargin < 2
%         force = 0;
%     end
%     
%     % if channels aren't specified
%     if nargin < 3
%         channels = 1:state.mcViewer.totalChannels;
%     end
% 
%     % data is already filtered:
%     if ~force && size(state.mcViewer.tsFilteredData, 2) >= fileCounter && ~isempty(state.mcViewer.tsFilteredData{1, fileCounter})
%         return
%     end
%     
%     % no data is loaded to filter
%     if isempty(state.mcViewer.tsData)
%         return
%     end
% 
% 
% 
%     %% Filter Data- first mcChannel Data
%     
%     % if any mcChannel is specified filter them all and continue to
%     % auxiliary channels.  mcChannel filtering is controlled separately
%     % from channelControl GUI elements.  ChannelConrtrol GUI elements only
%     % determine display filtering for these channels.
%     if any(channels < state.mcViewer.nChannels)
%         lowCutoff = 0;
%         highCutoff = 0;
%         if state.mcViewer.lowPass
%             lowCutoff = state.mcViewer.lowCutoff;
%         end
%         if state.mcViewer.highPass
%             highCutoff = state.mcViewer.highCutoff;
%         end
%  
%         filtData = mcFilter(state.mcViewer.tsData{1, fileCounter}(:,1:state.mcViewer.nChannels), lowCutoff, highCutoff, state.mcViewer.sampleRate);
%         state.mcViewer.analysis.tsSD(:,fileCounter) = std(filtData, 0, 1)';
%         state.mcViewer.analysis.tsMeans(:, fileCounter) = mean(filtData, 1)';
%         if state.mcViewer.blankTrain
%             pulseTrain = repmat(mcMakePulseTrain', 1, state.mcViewer.nChannels);
%             if isempty(pulseTrain)
%                 state.mcViewer.tsFilteredData{1, fileCounter}(:, 1:state.mcViewer.nChannels) = filtData;
%             else
%                 state.mcViewer.tsFilteredData{1, fileCounter}(:, 1:state.mcViewer.nChannels) = filtData .* pulseTrain; %FS MOD Kludge- need to add GUI controls for this...
%             end
%         else
%             state.mcViewer.tsFilteredData{1, fileCounter}(:, 1:state.mcViewer.nChannels) = filtData;
%         end
%     end
%     %% Filter Data- second, auxiliary channels
%     
%     % if aux channels don't exist, return
%     if ~state.mcViewer.totalChannels > state.mcViewer.nChannels
%         return
%     end
%     
%     for i = state.mcViewer.nChannels + 1 : state.mcViewer.totalChannels
%         if state.mcViewer.channel(i).ShowFilter
%             % filter the channel
%             state.mcViewer.tsFilteredData{1, fileCounter}(:, i) = ...
%                 mcFilter(state.mcViewer.tsData{1, fileCounter}(:, i), ...
%                 state.mcViewer.channel(i).LowPass, ...
%                 state.mcViewer.channel(i).HighPass ...
%                 );          
%         else
%             state.mcViewer.tsFilteredData{1, fileCounter}(:, i) = ...
%                             state.mcViewer.tsData{1, fileCounter}(:, i);
%         end
%     end
%      