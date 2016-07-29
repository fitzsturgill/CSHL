global gh state
format long;
h=waitbar(0, 'Setting Up Image Analysis Environment...', 'name', 'Image Processing Set-up');
waitbar(1/3, h,'Setting-up figures...');
gh.imageProcessingGUI = guihandles(imageProcessingGUI);
gh.imageParsingGUI = guihandles(imageParsingGUI);
gh.maxProjectionGUI = guihandles(maxProjectionGUI);
gh.montageGUI = guihandles(montageGUI);
gh.fileCounterGUI = guihandles(fileCounterGUI);
gh.overlayGUI = guihandles(overlayGUI);
gh.movieGUI = guihandles(movieGUI);
gh.averagingGUI = guihandles(averagingGUI);
gh.loadArrayGUI = guihandles(loadArrayGUI);
gh.genericLUTGUI = guihandles(genericLUTGUI);
gh.spineGUI = guihandles(spineGUI);
gh.spineDataGUI = guihandles(spineDataGUI);
gh.autotransformGUI = guihandles(autotransformGUI);
gh.excelLinkGUI = guihandles(excelLinkGUI);
gh.mathGUI = guihandles(mathGUI);
gh.roiAnalysis=guihandles(roiAnalysis);

waitbar(2/3, h,'Initializing...');
disp('Reading initialization file...');
openini('imageInit.ini');
waitbar(4/5, h,'Finishing Set-up...');
switchManualAuto;
makeSpineAnalysisFigures;
resetSpines;
resetDendrite;
state.imageProc.spine.imagehandle = image('CData', zeros(512,512), ...
	'CDataMapping', 'scaled', 'Parent', gh.spineGUI.initialaxis, 'ButtonDownFcn', 'spineImageOverFcn');
		state.imageProc.spine.mainImage = image('CData', zeros(512,512),...
			'CDataMapping', 'scaled','Parent', gh.spineGUI.mainAxes);


waitbar(1, h,'Done');
close(h);
