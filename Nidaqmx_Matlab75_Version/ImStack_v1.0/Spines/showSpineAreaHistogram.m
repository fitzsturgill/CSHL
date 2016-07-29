function showSpineAreaHistogram
global state gh

try
	figure, hist(state.imageProc.spine.globalSpineArea,state.imageProc.spine.histogramBins);
catch
	disp('Must add spine stats to list');
end
