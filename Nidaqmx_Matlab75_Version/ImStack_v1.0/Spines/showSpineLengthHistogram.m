function showSpineLengthHistogram
global state gh

try
	figure, hist(state.imageProc.spine.globalSpineLength,state.imageProc.spine.histogramBins);
catch
	disp('Must add spine stats to list');
end