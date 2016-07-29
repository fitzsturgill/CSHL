function showCumLengthHistogram(type)
global gh state

if nargin < 1
	figure('name', ['Cumulative Spine Length Histogram for ' state.imageProc.spineData.spineDataName], 'NumberTitle', 'off');
	plot(state.imageProc.spineData.HistogramData(:,1),makeCumFromRaw(state.imageProc.spineData.HistogramData(:,2)),'LineWidth',2)
	set(get(gca,'YLabel'),'String','Fraction of Spines with Length Less than X');
	set(get(gca,'XLabel'),'String','Spine Length (um)');
		set(gca, 'YLim', [0 1.4]);
	title( ['Cumulative Spine Length Histogram for ' state.imageProc.spineData.spineDataName]);

else
	if strcmp(type, 'all')
		figure('name', 'Cumulative Spine Length Histogram for All Files', 'NumberTitle', 'off');
		for i = 1:size(state.imageProc.spineData.allDataList,1)
			data(:,i) = state.imageProc.spineData.allDataList{i,3}(:,2);
		end
		
		plot(state.imageProc.spineData.HistogramData(:,1),makeCumFromRaw(data),'LineWidth',2);
		set(get(gca,'YLabel'),'String','Fraction of Spines with Length Less than X');
		set(get(gca,'XLabel'),'String','Spine Length (um)');
		title( ['Cumulative Spine Length Histogram for All Files']);
		legend(gca,state.imageProc.spineData.allDataList{:,1});
			set(gca, 'YLim', [0 1.4]);
	end
end
set(gcf, 'Color' , [1 1 1]);
copyToClip(gcf);