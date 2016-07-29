function saveDataInFigure(chan)
	if nargin<1
		chan=1;
	end
	global state

	axisPosition = [0 0 1 1];
	aspectRatio = state.internal.imageAspectRatioBias*state.acq.scanAmplitudeY/state.acq.scanAmplitudeX; 
		
	for channelCounter = chan
		low = getfield(state.internal, ['lowPixelValue' num2str(channelCounter)]);
		high = getfield(state.internal, ['highPixelValue' num2str(channelCounter)]);
		
		thisFigure = figure(...
			'Position', eval(['state.windowPositions.image' num2str(channelCounter) '_position']), ...
			'doublebuffer', 'on',...
			'KeyPressFcn', 'genericKeyPressFunction', ...
			'Tag',  ['GraphFigure' num2str(channelCounter)], ...
			'Name',  ['Saved Channel ' num2str(channelCounter) ' acq ' num2str(state.files.fileCounter-1)], ...
			'NumberTitle', 'off',  ...
			'MenuBar', 'figure', ...
			'CloseRequestFcn', 'set(gcf, ''visible'', ''off'')', ...
			'visible', 'on', ...
			'ColorMap', makeColorMap('graySat',8) ...
			);
		
		thisAxis = axes(...
			'YDir', 'Reverse', ...
			'NextPlot', 'add', ...
			'XLim', [0 state.acq.pixelsPerLine], ...
			'YLim', [0 state.acq.linesPerFrame], ...
			'CLim', [low high], ...
			'Parent', thisFigure, ...
			'YTickLabelMode', 'manual', ...
			'XTickLabelMode', 'manual', ...
			'Position', axisPosition,  ...
			'XTickLabel', [], ...
			'YTickLabel', [], ...
			'DataAspectRatioMode', 'manual',  ...
			'DataAspectRatio', [aspectRatio 1 1], ...
			'ButtonDownFcn', 'figureButtonOverCallback'...
			);
		
		thisImage = image(...
			state.acq.acquiredData{channelCounter}, ...
			'CDataMapping', 'Scaled', ...
			'Parent', thisAxis, ...
			'ButtonDownFcn', 'figureButtonOverCallback'...
			);
    end
    
    %FS Fitz MOD
%             %replicate blaster Locations on saved image given that
%        % ReferenceFigure is being used to select blaster positions
%         refHandles=get(state.internal.refAxis, 'children'); 
%         for i=1:length(refHandles)
%             child=get(refHandles(i));
%             fields=fieldnames(child);
%             if strcmp('text', child.Type)
%                 newTextHandles(i)=text();
%                 for j=1:length(fields)
%                     try
%                         set(newTextHandles(i), fields{j}, child.(fields{j}));
%                     end
%                     set(newTextHandles(i), 'Parent', thisAxis)
%                 end
%             elseif strcmp('line', child.Type)
%                 newLineHandles(i)=line(); 
%                 for j=1:length(fields)
%                     try
%                         set(newLineHandles(i), fields{j}, child.(fields{j}));
%                     end
%                     set(newLineHandles(i), 'Parent', thisAxis)                    
%                 end
%             end
%         end
% 	end
