function switchManualAuto
global gh state

if state.imageProc.spine.auto
	set(gh.spineGUI.autoSpineAnalysis, 'enable', 'on');
	set(gh.spineGUI.manualSpineAnalysis, 'enable', 'off');
	set(gh.spineGUI.text9, 'Visible', 'on');
	set(gh.spineGUI.text19, 'Visible', 'on');
	set(gh.spineGUI.text7, 'Visible', 'on');
	set(gh.spineGUI.text23, 'Visible', 'on');
	set(gh.spineGUI.text21, 'Visible', 'on');
	set(gh.spineGUI.text22, 'Visible', 'on');
	set(gh.spineGUI.frame5, 'Visible', 'on');
	set(gh.spineGUI.spineThickness, 'Visible', 'on');
	set(gh.spineGUI.spineThicknessSlider, 'Visible', 'on');
	set(gh.spineGUI.spineThreshold, 'Visible', 'on');
	set(gh.spineGUI.dendriteThickness, 'Visible', 'on');
	set(gh.spineGUI.dendriteThicknessSlider, 'Visible', 'on');
	set(gh.spineGUI.totalLength, 'Visible', 'on');
	set(gh.spineGUI.totalSpines, 'Visible', 'on');
	set(gh.spineGUI.totalDensity, 'Visible', 'on');
	set(gh.spineGUI.edit31, 'Visible', 'off');
	set(gh.spineGUI.edit32, 'Visible', 'off');
	set(gh.spineGUI.text33, 'Visible', 'off');
	set(gh.spineGUI.row, 'Visible', 'off');
	set(gh.spineGUI.column, 'Visible', 'off');
	set(gh.spineGUI.text39, 'Visible', 'off');
	set(gh.spineGUI.text40, 'Visible', 'off');
	set(gh.spineGUI.text41, 'Visible', 'off');
	set(gh.spineGUI.text42, 'Visible', 'off');
	set(gh.spineGUI.frame10, 'Visible', 'off');
	set(gh.spineGUI.pathName, 'Visible', 'off');
	set(gh.spineGUI.baseName, 'Visible', 'off');
	set(gh.spineGUI.savePreview, 'Visible', 'off');
	set(gh.spineGUI.fileSave, 'Visible', 'off');
	set(gh.spineGUI.fileSaveSlider, 'Visible', 'off');
	set(gh.spineGUI.frame12, 'Visible', 'off');
	set(gh.spineGUI.histogramBins, 'Position', [153.40000000000001   5  8.0   1.38461538461538]);
	set(gh.spineGUI.histogramBinsSlider, 'Position', [153.40000000000001   4.15  8.0 0.69230769230769]);
	set(gh.spineGUI.text18, 'Position', [136.6 5  15.600000000000001   1.153846153846154]);
	
else
	set(gh.spineGUI.autoSpineAnalysis, 'enable', 'off');
	set(gh.spineGUI.manualSpineAnalysis, 'enable', 'on');
	set(gh.spineGUI.text9, 'Visible', 'off');
	set(gh.spineGUI.text19, 'Visible', 'off');
	set(gh.spineGUI.text7, 'Visible', 'off');
	set(gh.spineGUI.edit31, 'Visible', 'on');
	set(gh.spineGUI.edit32, 'Visible', 'on');
	set(gh.spineGUI.text33, 'Visible', 'on');
	set(gh.spineGUI.text23, 'Visible', 'off');
	set(gh.spineGUI.text21, 'Visible', 'off');
	set(gh.spineGUI.text22, 'Visible', 'off');
	set(gh.spineGUI.frame5, 'Visible', 'off');
	set(gh.spineGUI.row, 'Visible', 'on');
	set(gh.spineGUI.column, 'Visible', 'on');
	set(gh.spineGUI.spineThreshold, 'Visible', 'off');
	set(gh.spineGUI.spineThickness, 'Visible', 'off');
	set(gh.spineGUI.spineThicknessSlider, 'Visible', 'off');
	set(gh.spineGUI.dendriteThickness, 'Visible', 'off');
	set(gh.spineGUI.dendriteThicknessSlider, 'Visible', 'off');
	set(gh.spineGUI.totalLength, 'Visible', 'off');
	set(gh.spineGUI.totalSpines, 'Visible', 'off');
	set(gh.spineGUI.totalDensity, 'Visible', 'off');
	set(gh.spineGUI.text39, 'Visible', 'on');
	set(gh.spineGUI.text40, 'Visible', 'on');
	set(gh.spineGUI.text41, 'Visible', 'on');
	set(gh.spineGUI.text42, 'Visible', 'on');
	set(gh.spineGUI.frame10, 'Visible', 'on');
	set(gh.spineGUI.pathName, 'Visible', 'on');
	set(gh.spineGUI.baseName, 'Visible', 'on');
	set(gh.spineGUI.savePreview, 'Visible', 'on');
	set(gh.spineGUI.fileSave, 'Visible', 'on');
	set(gh.spineGUI.fileSaveSlider, 'Visible', 'on');
	set(gh.spineGUI.frame12, 'Visible', 'on');
	set(gh.spineGUI.histogramBins, 'Position', [154.6   1.423  8.0  1.38461538461538]);
	set(gh.spineGUI.histogramBinsSlider, 'Position', [154.6   .677   8.0   0.69230769230769]);
	set(gh.spineGUI.text18, 'Position', [137  1.677   15.6 1.153846153846154]);
end
