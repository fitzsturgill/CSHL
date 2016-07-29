function loadImageFromArray(imageIn)
global gh state

if state.imageProc.internal.imageCounter < 1
	state.imageProc.internal.imageCounter = 1;
end
state.imageProc.internal.loadArrayCounter = state.imageProc.internal.loadArrayCounter+1;
pname = get(gh.fileCounterGUI.pathName, 'String');
filename = [pname get(gh.fileCounterGUI.baseName, 'String')];
filename = [getbaseName(filename,pname) num2str(state.imageProc.internal.loadArrayCounter)];
try	
	cd(pname);
catch
end

resetSlidersToMax;

eval(['state.imageProc.cell.currentImage{state.imageProc.internal.imageCounter} = ' imageIn ';']);
updateGUIByGlobalCell('state.imageProc.currentImage', state.imageProc.internal.imageCounter);

% get the size of the image (total number of frames)
state.imageProc.cell.numberOfFrames{state.imageProc.internal.imageCounter} = ...
	size(state.imageProc.cell.currentImage{state.imageProc.internal.imageCounter},3);
updateGUIByGlobalCell('state.imageProc.numberOfFrames', state.imageProc.internal.imageCounter);

% get the pixels per line (X) of the image 
state.imageProc.parsing.cell.pixelsPerLine{state.imageProc.internal.imageCounter} = ...
	size(state.imageProc.cell.currentImage{state.imageProc.internal.imageCounter}, 2);
updateGUIByGlobalCell('state.imageProc.parsing.pixelsPerLine', state.imageProc.internal.imageCounter);

% get the lines per frame (Y)of the image 
state.imageProc.parsing.cell.linesPerFrame{state.imageProc.internal.imageCounter} ...
	= size(state.imageProc.cell.currentImage{state.imageProc.internal.imageCounter},1);
updateGUIByGlobalCell('state.imageProc.parsing.linesPerFrame', state.imageProc.internal.imageCounter);


% update montage end
state.imageProc.cell.montageEnd{state.imageProc.internal.imageCounter} = state.imageProc.cell.numberOfFrames ...
	{state.imageProc.internal.imageCounter};
updateGUIByGlobalCell('state.imageProc.montageEnd', state.imageProc.internal.imageCounter);

% update max end
state.imageProc.cell.maxEnd{state.imageProc.internal.imageCounter} = state.imageProc.cell.numberOfFrames ...
	{state.imageProc.internal.imageCounter};
updateGUIByGlobalCell('state.imageProc.maxEnd', state.imageProc.internal.imageCounter);

% update movie end
state.imageProc.cell.movieEnd{state.imageProc.internal.imageCounter} = state.imageProc.cell.numberOfFrames ...
	{state.imageProc.internal.imageCounter};
updateGUIByGlobalCell('state.imageProc.movieEnd', state.imageProc.internal.imageCounter);

% update the current frame in the image
state.imageProc.cell.currentFrame{state.imageProc.internal.imageCounter} = 1;
updateGUIByGlobalCell('state.imageProc.currentFrame', state.imageProc.internal.imageCounter);

% update file name in popupmenu
% Update the file counters also if they are to be updated
% value of popupmenu is the counter....all logic based on this.
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
	updatePopUpMenuStringByValue(gh.imageProcessingGUI.fileName, state.imageProc.internal.imageCounter, filename);
else
	updateNextPopUpWithString(gh.imageProcessingGUI.fileName, filename);
end

% Set lookup table to scale for image

state.imageProc.cell.lowPixelValue{state.imageProc.internal.imageCounter} = state.imageProc.lowPixelValue;
updateGUIByGlobalCell('state.imageProc.lowPixelValue', state.imageProc.internal.imageCounter);

state.imageProc.cell.highPixelValue{state.imageProc.internal.imageCounter} = state.imageProc.highPixelValue;
updateGUIByGlobalCell('state.imageProc.highPixelValue', state.imageProc.internal.imageCounter);


% Parse the image.
% Insert Bernardo's Parsing Function here

state.imageProc.cell.currentFrame{state.imageProc.internal.imageCounter} = 1;
state.imageProc.cell.montageStart{state.imageProc.internal.imageCounter} = 1;
state.imageProc.cell.maxStart{state.imageProc.internal.imageCounter} = 1;
state.imageProc.cell.movieStart{state.imageProc.internal.imageCounter} = 1;
state.imageProc.parsing.cell.numberOfChannels{state.imageProc.internal.imageCounter} = 1;
state.imageProc.parsing.cell.numberOfFrames{state.imageProc.internal.imageCounter} = 1;
state.imageProc.parsing.cell.scanRotation{state.imageProc.internal.imageCounter} = 0;
state.imageProc.parsing.cell.averaged{state.imageProc.internal.imageCounter} = 0;
state.imageProc.parsing.cell.numberOfZSlices{state.imageProc.internal.imageCounter} = state.imageProc.numberOfFrames;
updateGUIByGlobalCell('state.imageProc.currentFrame', state.imageProc.internal.imageCounter);
updateGUIByGlobalCell('state.imageProc.montageStart', state.imageProc.internal.imageCounter);
updateGUIByGlobalCell('state.imageProc.maxStart', state.imageProc.internal.imageCounter);
updateGUIByGlobalCell('state.imageProc.movieStart', state.imageProc.internal.imageCounter);
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

% update filenameGUI popupmenu
state.imageProc.fileNameGUI = state.imageProc.internal.imageCounter;
updateGUIByGlobal('state.imageProc.fileNameGUI');

%Update other Popupmenus with file Name Values.
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
		'''imageProcKeyPressF'',''Tag'', filename, ''Name'', '...
		'filename,''MenuBar'', ''figure'',''NumberTitle'', ''off'','...
		'''Position'', positionF,''ButtonDownFcn'',''figureButtonDownFcn'', ''CloseRequestFcn'', ''figureCloseFunction'',''NextPlot'', ''Add'');' ]);

colormap(makeColorMap('gray',8));

eval(['state.imageProc.internal.axis{state.imageProc.internal.imageCounter} = axes(''YDir'', ''Reverse'',''NextPlot'', ''Add'', ''XLim'', '...
		'[0 state.imageProc.parsing.cell.pixelsPerLine{state.imageProc.internal.imageCounter} ],''YLim'', [0 ' ...
		'state.imageProc.parsing.cell.linesPerFrame{state.imageProc.internal.imageCounter}],''CLim'', '...
		'[state.imageProc.cell.lowPixelValue{state.imageProc.internal.imageCounter} state.imageProc.cell.highPixelValue{state.imageProc.internal.imageCounter}],'...
		'''Parent'', state.imageProc.internal.Figure{state.imageProc.internal.imageCounter}, ''DataAspectRatioMode'', ''manual'', ''YTickLabelMode'', ''manual'',' ...
		' ''XTickLabelMode'', ''manual'', ''XTickLabel'', [] , ''YTickLabel'', [], ''Position'', [0 0 1 1]);']);

eval(['state.imageProc.internal.imagehandle{state.imageProc.internal.imageCounter} = image(''CData'', state.imageProc.cell.currentImage' ...
		'{state.imageProc.internal.imageCounter}(:,:,state.imageProc.cell.currentFrame{state.imageProc.internal.imageCounter}),''CDataMapping'', ''Scaled'', ''Parent'', '...
		'state.imageProc.internal.axis{state.imageProc.internal.imageCounter}, ''ButtonDownFcn'', ''figureButtonOverCallbackImagePr'');']);

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
