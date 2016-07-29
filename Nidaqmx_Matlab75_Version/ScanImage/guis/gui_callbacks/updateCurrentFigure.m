function updateCurrentFigure
global state

% updateCurrentFigure.m****
% Function that will grab the current figure positions and write them to the structure state.
%
% Written By: Thomas Pologruto
% Cold Spring Harbor Labs
% March 2, 2001

imageCounter = 0;

try
	
	for channelCounter = 1:state.init.maximumNumberOfInputChannels
		channelOn = eval(['state.acq.imagingChannel' num2str(channelCounter)]);
		if channelOn == 1	
			imageCounter = imageCounter + 1;
			eval(['position{imageCounter} = get(state.internal.GraphFigure' num2str(imageCounter) ', ''Position'');']);
			eval(['state.internal.figurePositionX' num2str(imageCounter) ' = position{imageCounter}(1,1);']);
			eval(['state.internal.figurePositionY' num2str(imageCounter) '= position{imageCounter}(1,2);']);
			eval(['state.internal.figureWidth' num2str(imageCounter) ' = position{imageCounter}(1,3);']);
			eval(['state.internal.figureHeight' num2str(imageCounter) ' = position{imageCounter}(1,4);']);
			end
	end

	imageCounter = 0;

	for channelCounter = 1:state.init.maximumNumberOfInputChannels
		channelOn = eval(['state.acq.maxImage' num2str(channelCounter)]);
		if channelOn == 1	
			imageCounter = imageCounter + 1;
			eval(['position{imageCounter} = get(state.internal.MaxFigure' num2str(imageCounter) ', ''Position'');']);
			eval(['state.internal.maxfigurePositionX' num2str(imageCounter) ' = position{imageCounter}(1,1);']);
			eval(['state.internal.maxfigurePositionY' num2str(imageCounter) '= position{imageCounter}(1,2);']);
			eval(['state.internal.maxfigureWidth' num2str(imageCounter) ' = position{imageCounter}(1,3);']);
			eval(['state.internal.maxfigureHeight' num2str(imageCounter) ' = position{imageCounter}(1,4);']);
		end
	end
catch
	disp(['updateCurrentFigure: ' lasterr]);
end
