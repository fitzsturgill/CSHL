function image2bw(array, spineThreshold)
global gh state

array = array(:,:,state.imageProc.spine.currentSpineFrame) < spineThreshold;
state.imageProc.spine.bw = ~array;

columns = size(array,2);
rows = size(array,1);
set(gh.spineGUI.mainAxes, 'YLim', [1 rows], 'XLim', [1 columns], 'Ydir', 'reverse', 'Clim', [0 1]);
set(gh.spineGUI.figure1, 'Colormap', makeColormap('gray',8));
state.imageProc.spine.bwimagehandle = image('CData', state.imageProc.spine.bw, 'CDataMapping', 'scaled','Parent', gh.spineGUI.mainAxes);




