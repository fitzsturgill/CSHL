function selectSpinesManually_HW
global gh state bit
% bit is 0 if it is a normal point and 1 if it si an extension.
bit = [];
peaks=[];

figure(gh.spineGUI.mainFigure);
data=smooth2(improfile', 3);

offset=min(data);

waveo('FWHMdata', data);
[pd, py]=findpeaks_small(data, 1, offset, 0.5);
waveo('FWHMx', pd-1);
waveo('FWHMy', py);
hw=pd(2)-pd(1);
	
disp(num2str(hw));

state.imageProc.spine.numberOfSpines = state.imageProc.spine.numberOfSpines + 1;
updateGUIByGlobal('state.imageProc.spine.numberOfSpines');


state.imageProc.spine.spineLengths(end+1) = hw;

%saveSpineDataToExcel;