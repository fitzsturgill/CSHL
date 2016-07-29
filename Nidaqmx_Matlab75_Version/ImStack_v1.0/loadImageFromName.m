function loadImageFromName(filename,pname)
global gh state


if state.imageProc.internal.imageCounter < 1
	state.imageProc.internal.imageCounter = 1;
end

try	
	cd(pname);
catch
end

resetSlidersToMax;

% set the first current image array to the first image opened.
[path,name,ext] = fileparts(filename);
switch ext
case '.tif'
	state.imageProc.cell.currentImage{state.imageProc.internal.imageCounter} = opentif(filename);
case '.CFD'
	state.imageProc.cell.currentImage{state.imageProc.internal.imageCounter} = openACFDSpine(filename, state.imageProc.cfd.numberofChannels);
end

if ~isnumeric(state.imageProc.cell.currentImage{state.imageProc.internal.imageCounter})
	return
else
end

updateGUIByGlobalCell('state.imageProc.currentImage', state.imageProc.internal.imageCounter);

% get the size of the image (total number of frames)
if state.imageProc.colorMap == 0
	state.imageProc.cell.numberOfFrames{state.imageProc.internal.imageCounter} = ...
		size(state.imageProc.cell.currentImage{state.imageProc.internal.imageCounter},3);
	updateGUIByGlobalCell('state.imageProc.numberOfFrames', state.imageProc.internal.imageCounter);
else
	state.imageProc.cell.numberOfFrames{state.imageProc.internal.imageCounter} = 1;
	updateGUIByGlobalCell('state.imageProc.numberOfFrames', state.imageProc.internal.imageCounter);
end

% get the pixels per line (X) of the image 
state.imageProc.parsing.cell.pixelsPerLine{state.imageProc.internal.imageCounter} = ...
	size(state.imageProc.cell.currentImage{state.imageProc.internal.imageCounter}, 2);
updateGUIByGlobalCell('state.imageProc.parsing.pixelsPerLine', state.imageProc.internal.imageCounter);

% get the lines per frame (Y)of the image 
state.imageProc.parsing.cell.linesPerFrame{state.imageProc.internal.imageCounter} ...
	= size(state.imageProc.cell.currentImage{state.imageProc.internal.imageCounter},1);
updateGUIByGlobalCell('state.imageProc.parsing.linesPerFrame', state.imageProc.internal.imageCounter);

% update movie end
state.imageProc.cell.movieEnd{state.imageProc.internal.imageCounter} = state.imageProc.cell.numberOfFrames ...
	{state.imageProc.internal.imageCounter};
updateGUIByGlobalCell('state.imageProc.movieEnd', state.imageProc.internal.imageCounter);

% update movie start
state.imageProc.cell.movieStart{state.imageProc.internal.imageCounter} = 1;
updateGUIByGlobalCell('state.imageProc.movieStart', state.imageProc.internal.imageCounter);

% update montage end
state.imageProc.cell.montageEnd{state.imageProc.internal.imageCounter} = state.imageProc.cell.numberOfFrames ...
	{state.imageProc.internal.imageCounter};
updateGUIByGlobalCell('state.imageProc.montageEnd', state.imageProc.internal.imageCounter);

% update montage start
state.imageProc.cell.montageStart{state.imageProc.internal.imageCounter} = 1;
updateGUIByGlobalCell('state.imageProc.montageStart', state.imageProc.internal.imageCounter);

% update max end
state.imageProc.cell.maxEnd{state.imageProc.internal.imageCounter} = state.imageProc.cell.numberOfFrames ...
	{state.imageProc.internal.imageCounter};
updateGUIByGlobalCell('state.imageProc.maxEnd', state.imageProc.internal.imageCounter);

% update max start
state.imageProc.cell.maxStart{state.imageProc.internal.imageCounter} = 1;
updateGUIByGlobalCell('state.imageProc.maxStart', state.imageProc.internal.imageCounter);

% update the current frame in the image
state.imageProc.cell.currentFrame{state.imageProc.internal.imageCounter} = 1;
updateGUIByGlobalCell('state.imageProc.currentFrame', state.imageProc.internal.imageCounter);

% update file name in popupmenu
% Update the file counters also if they are to be updated
% value of popupmenu is the counter....all logic based on this.
str = get(gh.imageProcessingGUI.fileName, 'String');
if ischar(str)
	if strcmp(str, filename)
		filename = [filename num2str(state.imageProc.internal.nameCounter)];
	end
elseif iscell(str)
	for i = 1:length(str)
		if strcmp(str{i}, filename)
			filename = [filename num2str(state.imageProc.internal.nameCounter)];
		end
	end
end

state.imageProc.cell.fileName{state.imageProc.internal.imageCounter} = filename;
updateGUIByGlobalCell('state.imageProc.fileName', state.imageProc.internal.imageCounter);
state.imageProc.parsing.cell.header{state.imageProc.internal.imageCounter} = readImageHeaderTif(filename);
updateGUIByGlobalCell('state.imageProc.parsing.header', state.imageProc.internal.imageCounter);

if state.imageProc.internal.imageCounter == 1 & state.imageProc.auto.loadAsTransform == 0 & state.imageProc.internal.updateFileInfo == 1% First image loaded
	
	state.imageProc.pathName = pname;
	state.imageProc.baseName= getbasename(filename, pname);
	updateGUIByGlobal('state.imageProc.baseName');
	updateGUIByGlobal('state.imageProc.pathName');
	
elseif state.imageProc.internal.imageCounter > 1 & state.imageProc.internal.updateFileInfo == 1
	
	state.imageProc.pathName = pname;
	state.imageProc.baseName= getbasename(filename, pname);
	updateGUIByGlobal('state.imageProc.baseName');
	updateGUIByGlobal('state.imageProc.pathName');
	
elseif state.imageProc.internal.imageCounter > 1 & state.imageProc.internal.updateFileInfo == 0
end

if state.imageProc.internal.imageCounter == 1
	updatePopUpMenuStringByValue(gh.imageProcessingGUI.fileName, state.imageProc.internal.imageCounter, filename)
else
	updateNextPopUpWithString(gh.imageProcessingGUI.fileName, filename);
end

% Set lookup table to scale for image
try
	
	if state.imageProc.updateClim == 1 
		
		state.imageProc.cell.lowPixelValue{state.imageProc.internal.imageCounter} = min(min(min(state.imageProc.cell.currentImage ... 
			{state.imageProc.internal.imageCounter})));
		updateGUIByGlobalCell('state.imageProc.lowPixelValue', state.imageProc.internal.imageCounter);
		
		state.imageProc.cell.highPixelValue{state.imageProc.internal.imageCounter} = (.01*state.imageProc.LUTpercentmax)*double(max(max(max(state.imageProc.cell.currentImage ...
			{state.imageProc.internal.imageCounter}))));
		updateGUIByGlobalCell('state.imageProc.highPixelValue', state.imageProc.internal.imageCounter);
		
	else 
		
		state.imageProc.cell.lowPixelValue{state.imageProc.internal.imageCounter} = state.imageProc.lowPixelValue;
		updateGUIByGlobalCell('state.imageProc.lowPixelValue', state.imageProc.internal.imageCounter);
		
		state.imageProc.cell.highPixelValue{state.imageProc.internal.imageCounter} = state.imageProc.highPixelValue;
		updateGUIByGlobalCell('state.imageProc.highPixelValue', state.imageProc.internal.imageCounter);
	end
	
	if state.imageProc.lowPixelValue >= state.imageProc.highPixelValue 
		state.imageProc.cell.highPixelValue{state.imageProc.internal.imageCounter} = double(state.imageProc.highPixelValue) + 1;
		updateGUIByGlobalCell('state.imageProc.highPixelValue', state.imageProc.internal.imageCounter);
	end
	
catch
	state.imageProc.cell.lowPixelValue{state.imageProc.internal.imageCounter} = 0;
	updateGUIByGlobalCell('state.imageProc.lowPixelValue', state.imageProc.internal.imageCounter);
	
	state.imageProc.cell.highPixelValue{state.imageProc.internal.imageCounter} = 100;
	updateGUIByGlobalCell('state.imageProc.highPixelValue', state.imageProc.internal.imageCounter);
end

state.imageProc.cell.currentFrame{state.imageProc.internal.imageCounter} = 1;
state.imageProc.cell.montageStart{state.imageProc.internal.imageCounter} = 1;
state.imageProc.cell.maxStart{state.imageProc.internal.imageCounter} = 1;

% Parse the image.
% These dont need parsing
% Parse the header if possible

if strcmp('No Image Description Found.', state.imageProc.parsing.header) %No header
	state.imageProc.parsing.cell.numberOfChannels{state.imageProc.internal.imageCounter} = 1;
	state.imageProc.parsing.cell.numberOfFrames{state.imageProc.internal.imageCounter} = 1; 
	state.imageProc.parsing.cell.scanRotation{state.imageProc.internal.imageCounter} = 0;
	state.imageProc.parsing.cell.averaged{state.imageProc.internal.imageCounter} = 0;
	state.imageProc.parsing.cell.numberOfZSlices{state.imageProc.internal.imageCounter} = 1;
else
	try % Try reading ScanImage Header
		state.imageProc.parsing.cell.numberOfChannels{state.imageProc.internal.imageCounter} = ...
			valueFromHeaderString('state.acq.numberOfChannelsSave', state.imageProc.parsing.header);
		state.imageProc.parsing.cell.numberOfFrames{state.imageProc.internal.imageCounter} = ...
			valueFromHeaderString('state.acq.numberOfFrames', state.imageProc.parsing.header);
		state.imageProc.parsing.cell.scanRotation{state.imageProc.internal.imageCounter} = ...
			valueFromHeaderString('state.acq.scanRotation', state.imageProc.parsing.header);
		state.imageProc.parsing.cell.averaged{state.imageProc.internal.imageCounter} = ...
			valueFromHeaderString('state.acq.averaging', state.imageProc.parsing.header);
		state.imageProc.parsing.cell.numberOfZSlices{state.imageProc.internal.imageCounter} = ...
			valueFromHeaderString('state.acq.numberOfZSlices', state.imageProc.parsing.header);
	catch
		try % In case you read another type of header (Flouview)
			numberOfChannels = 0;
			for i = 1:3
				try
					channelOn = valueFromHeaderString(['Channel ' num2str(i)], state.imageProc.parsing.header);
					numberofChannels = numberofChannels + 1;
				end
			end
			state.imageProc.parsing.cell.numberOfChannels{state.imageProc.internal.imageCounter} = numberOfChannels;
			average = valueFromHeaderString('Frame Filter', state.imageProc.parsing.header);
			
			if strcmp(' frame Kalman Filte', average)
				state.imageProc.parsing.cell.averaged{state.imageProc.internal.imageCounter} = 1;
				state.imageProc.parsing.cell.numberOfFrames{state.imageProc.internal.imageCounter} = 1;
				state.imageProc.parsing.cell.numberOfZSlices{state.imageProc.internal.imageCounter} = ...
					state.imageProc.cell.numberOfFrames{state.imageProc.internal.imageCounter} ;
				state.imageProc.parsing.cell.scanRotation{state.imageProc.internal.imageCounter} = 0;
			else
				state.imageProc.parsing.cell.averaged{state.imageProc.internal.imageCounter} = 0;
				state.imageProc.parsing.cell.numberOfFrames{state.imageProc.internal.imageCounter} = 1;
				state.imageProc.parsing.cell.numberOfZSlices{state.imageProc.internal.imageCounter} = ...
					state.imageProc.cell.numberOfFrames{state.imageProc.internal.imageCounter} ;
				state.imageProc.parsing.cell.scanRotation{state.imageProc.internal.imageCounter} = 0;
			end
		catch
			state.imageProc.parsing.cell.numberOfChannels{state.imageProc.internal.imageCounter} = 1;
			state.imageProc.parsing.cell.numberOfFrames{state.imageProc.internal.imageCounter} = 1; 
			state.imageProc.parsing.cell.scanRotation{state.imageProc.internal.imageCounter} = 0;
			state.imageProc.parsing.cell.averaged{state.imageProc.internal.imageCounter} = 0;
			state.imageProc.parsing.cell.numberOfZSlices{state.imageProc.internal.imageCounter} = 1;
		end
	end
end


updateGUIByGlobalCell('state.imageProc.currentFrame', state.imageProc.internal.imageCounter);
updateGUIByGlobalCell('state.imageProc.montageStart', state.imageProc.internal.imageCounter);
updateGUIByGlobalCell('state.imageProc.maxStart', state.imageProc.internal.imageCounter);
updateGUIByGlobalCell('state.imageProc.parsing.numberOfChannels', state.imageProc.internal.imageCounter);
updateGUIByGlobalCell('state.imageProc.parsing.numberOfFrames', state.imageProc.internal.imageCounter);
updateGUIByGlobalCell('state.imageProc.parsing.scanRotation', state.imageProc.internal.imageCounter);
updateGUIByGlobalCell('state.imageProc.parsing.averaged', state.imageProc.internal.imageCounter);
updateGUIByGlobalCell('state.imageProc.parsing.numberOfZSlices', state.imageProc.internal.imageCounter);


%Update other Popupmenus with file Name bars.
set(gh.maxProjectionGUI.fileName, 'String', get(gh.imageProcessingGUI.fileName, 'String'));
set(gh.montageGUI.fileName, 'String', get(gh.imageProcessingGUI.fileName, 'String'));
set(gh.movieGUI.fileName, 'String', get(gh.imageProcessingGUI.fileName, 'String'));
set(gh.overlayGUI.fileName1, 'String', get(gh.imageProcessingGUI.fileName, 'String')); 
set(gh.overlayGUI.fileName2, 'String', get(gh.imageProcessingGUI.fileName, 'String')); 
set(gh.overlayGUI.fileName3, 'String', get(gh.imageProcessingGUI.fileName, 'String')); 
set(gh.averagingGUI.fileName, 'String', get(gh.imageProcessingGUI.fileName, 'String')); 
set(gh.mathGUI.fileName1, 'String', get(gh.imageProcessingGUI.fileName, 'String')); 
set(gh.mathGUI.fileName2, 'String', get(gh.imageProcessingGUI.fileName, 'String')); 
set(gh.roiAnalysis.roiPopupmenu, 'String', get(gh.imageProcessingGUI.fileName, 'String')); 

%Update other Popupmenus with file Name Values.
% update filenameGUI popupmenu
state.imageProc.fileNameGUI = state.imageProc.internal.imageCounter;
updateGUIByGlobal('state.imageProc.fileNameGUI');

set(gh.maxProjectionGUI.fileName, 'Value', get(gh.imageProcessingGUI.fileName, 'Value'));
set(gh.montageGUI.fileName, 'Value', get(gh.imageProcessingGUI.fileName, 'Value'));
set(gh.movieGUI.fileName, 'Value', get(gh.imageProcessingGUI.fileName, 'Value'));
set(gh.overlayGUI.fileName1, 'Value', get(gh.imageProcessingGUI.fileName, 'Value')); 
set(gh.overlayGUI.fileName2, 'Value', get(gh.imageProcessingGUI.fileName, 'Value')); 
set(gh.overlayGUI.fileName3, 'Value', get(gh.imageProcessingGUI.fileName, 'Value')); 
set(gh.averagingGUI.fileName, 'Value', get(gh.imageProcessingGUI.fileName, 'Value')); 
set(gh.mathGUI.fileName1, 'Value', get(gh.imageProcessingGUI.fileName, 'Value')); 
set(gh.mathGUI.fileName2, 'Value', get(gh.imageProcessingGUI.fileName, 'Value')); 
set(gh.roiAnalysis.roiPopupmenu, 'Value', get(gh.imageProcessingGUI.fileName, 'Value')); 

%Update Mean for first frame
state.imageProc.internal.meanIntensity = mean2(state.imageProc.cell.currentImage{state.imageProc.internal.imageCounter}(:,:,state.imageProc.currentFrame));
updateGUIByGlobal('state.imageProc.internal.meanIntensity');
state.imageProc.internal.sumIntensity = sum(sum(state.imageProc.cell.currentImage{state.imageProc.internal.imageCounter}(:,:,state.imageProc.currentFrame)));
updateGUIByGlobal('state.imageProc.internal.sumIntensity');

% Select mode for image manipulations
% Creat figures, axes, and imagehandles.
if state.imageProc.parsing.pixelsPerLine > 512 | state.imageProc.parsing.linesPerFrame > 512
	positionF=[50 50  512 512];
else
	positionF=[50 50  1.05*state.imageProc.parsing.pixelsPerLine 1.05*state.imageProc.parsing.linesPerFrame];
end
eval(['state.imageProc.internal.Figure{state.imageProc.internal.imageCounter} = figure(''doublebuffer'', ''on'', ''KeyPressFcn'','...
		'''imageProcKeyPressF'',''Tag'', state.imageProc.cell.fileName{state.imageProc.internal.imageCounter}, ''Name'', '...
		'state.imageProc.cell.fileName{state.imageProc.internal.imageCounter},''MenuBar'', ''figure'',''NumberTitle'', ''off'','...
		'''Position'', positionF,''ButtonDownFcn'',''figureButtonDownFcn'', ''CloseRequestFcn'', ''figureCloseFunction'',''NextPlot'', ''Add'');' ]);

colormap(makeColorMap('gray',8));

if state.imageProc.lowPixelValue >= state.imageProc.highPixelValue 
	state.imageProc.cell.highPixelValue{state.imageProc.internal.imageCounter} = double(state.imageProc.lowPixelValue) + 1;
	updateGUIByGlobalCell('state.imageProc.highPixelValue', state.imageProc.internal.imageCounter);
end

eval(['state.imageProc.internal.axis{state.imageProc.internal.imageCounter} = axes(''YDir'', ''Reverse'',''NextPlot'', ''Add'', ''XLim'', '...
		'[0 state.imageProc.parsing.cell.pixelsPerLine{state.imageProc.internal.imageCounter} ],''YLim'', [0 ' ...
		'state.imageProc.parsing.cell.linesPerFrame{state.imageProc.internal.imageCounter}],''CLim'', '...
		'[state.imageProc.cell.lowPixelValue{state.imageProc.internal.imageCounter} state.imageProc.cell.highPixelValue{state.imageProc.internal.imageCounter}],'...
		'''Parent'', state.imageProc.internal.Figure{state.imageProc.internal.imageCounter}, ''DataAspectRatioMode'', ''manual'', ''YTickLabelMode'', ''manual'',' ...
		' ''XTickLabelMode'', ''manual'', ''XTickLabel'', [] , ''YTickLabel'', [], ''Position'', [0 0 1 1]);']);

if state.imageProc.colorMap == 0 % BW image
	eval(['state.imageProc.internal.imagehandle{state.imageProc.internal.imageCounter} = image(''CData'', state.imageProc.cell.currentImage' ...
			'{state.imageProc.internal.imageCounter}(:,:,state.imageProc.cell.currentFrame{state.imageProc.internal.imageCounter}),''CDataMapping'', ''Scaled'', ''Parent'', '...
			'state.imageProc.internal.axis{state.imageProc.internal.imageCounter}, ''ButtonDownFcn'', ''figureButtonOverCallbackImagePr'');']);
else
	eval(['state.imageProc.internal.imagehandle{state.imageProc.internal.imageCounter} = image(''CData'', state.imageProc.cell.currentImage' ...
			'{state.imageProc.internal.imageCounter},''CDataMapping'', ''Scaled'', ''Parent'', '...
			'state.imageProc.internal.axis{state.imageProc.internal.imageCounter}, ''ButtonDownFcn'', ''figureButtonOverCallbackImagePr'');']);
end

% Make UI Context Menus
cmenuHandle=createUiContextMenu('Save ROI', 'recordCurrentROI', 'Input ROI', 'inputCurrentROI', ...
	'Load Saved ROI','loadCurrentROI','Zoom To Saved ROI','zoomCurrentROI', 'Load and Process ROI',...
    'loadAndProcessCurROI','Histogram','histAndNoiseAnalysis','Time Series Mean','computeTimeSeriesMeanV',...
	'Pixel-by-Pixel Baseline', 'histAndNoiseAnalysisPixel', 'Threshold Pixel-byPixel','compareImgePixels','Copy','copyToClip(gco)',...
	'Average','averageAccrossZ','Project','projectAccrossZ','Median Filter', 'applyMedianFilter');
addUimenuToHandle(state.imageProc.internal.imagehandle{state.imageProc.internal.imageCounter}, cmenuHandle);

updateSliderMaxMin;

% increment the image counter
state.imageProc.internal.imageCounter = state.imageProc.internal.imageCounter + 1;
