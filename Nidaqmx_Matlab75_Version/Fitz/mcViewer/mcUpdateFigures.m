function mcUpdateFigures
    try
        global state

%         if state.phys.mcAcq.showThresh
%             show = 'on';
%         else
%             show = 'off';
%         end
        
        validChannels=find(mcChannelFieldToVector('Show')); %channels to be displayed in mcFigure
        state.phys.mcAcq.displayedChannels = validChannels;        
        nValidChannels = length(validChannels);
        for i = 1:nValidChannels
            channel = validChannels(i);
            set(state.phys.mcAcq.lines(1, i),...
                'Ydata', state.phys.mcAcq.displayData(:, channel),...
                'Xdata', state.phys.mcAcq.displayXData...
                );
        end
    catch
        lasterr
        disp('error in mcUpdateFigures');
    end
    
