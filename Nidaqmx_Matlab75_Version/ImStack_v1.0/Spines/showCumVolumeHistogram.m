function showCumVolumeHistogram(type)
global gh state

if nargin < 1
	figure('name', ['Cumulative Spine Volume Histogram for ' state.imageProc.spineData.spineDataName], 'NumberTitle', 'off');
	plot(state.imageProc.spineData.volumeHist(:,1),makeCumFromRaw(state.imageProc.spineData.volumeHist(:,2)),'LineWidth',2)
	set(get(gca,'YLabel'),'String','Fraction of Spines with Volume Less than X');
	set(get(gca,'XLabel'),'String','Spine Volume (fL)');
	set(gca, 'YLim', [0 1.4]);
	title( ['Cumulative Spine Volume Histogram for ' state.imageProc.spineData.spineDataName]);

else
	if strcmp(type, 'all')
		figure('name', 'Cumulative Spine Volume Histogram for All Files', 'NumberTitle', 'off');
		for i = 1:size(state.imageProc.spineData.allDataList,1)
			data(:,i) = state.imageProc.spineData.allDataList{i,4}(:,2);
		end
		
		plot(state.imageProc.spineData.volumeHist(:,1),makeCumFromRaw(data),'LineWidth',2);
		set(get(gca,'YLabel'),'String','Fraction of Spines with Volume Less than X');
		set(get(gca,'XLabel'),'String','Spine Volume (fL)');
		set(gca, 'YLim', [0 1.4]);
		title( ['Cumulative Spine Volume Histogram for All Files']);
		legend(gca,state.imageProc.spineData.allDataList{:,1});
	end
end
set(gcf, 'Color' , [1 1 1]);
copyToClip(gcf);