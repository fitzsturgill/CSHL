function showLengthHistogram(type)
global gh state

if nargin < 1
	figure('name', ['Spine Length Histogram for ' state.imageProc.spineData.spineDataName], 'NumberTitle', 'off');
	n = bar(state.imageProc.spineData.HistogramData(:,1), state.imageProc.spineData.HistogramData(:,2));
	set(get(gca,'YLabel'),'String','Number of Spines');
	set(get(gca,'XLabel'),'String','Spine Length (um)');
	title( ['Spine Length Histogram for ' state.imageProc.spineData.spineDataName]);
	h1 = gca;
	h2 = axes('Position',get(h1,'Position'));
	plot(state.imageProc.spineData.HistogramData(:,1),state.imageProc.spineData.HistogramData(:,2),'LineWidth',2)
	set(h2,'YAxisLocation','left','Color','none','XTickLabel',[],'YTickLabel', [])
	set(h2,'XLim',get(h1,'XLim'),'Layer','top')
	
else
	if strcmp(type, 'all')
		figure('name', 'Spine Length Histogram for All Files', 'NumberTitle', 'off');
		for i = 1:size(state.imageProc.spineData.allDataList,1)
			data(:,i) = state.imageProc.spineData.allDataList{i,3}(:,2);
		end
		n = bar(state.imageProc.spineData.HistogramData(:,1),data, 'group');
		set(get(gca,'YLabel'),'String','Number of Spines');
		set(get(gca,'XLabel'),'String','Spine Length (um)');
		title( ['Spine Length Histogram for All Files']);
		h1 = gca;
		h2 = axes('Position',get(h1,'Position'));
		plot(state.imageProc.spineData.HistogramData(:,1),data,'LineWidth',2)
		set(h2,'YAxisLocation','left','Color','none','XTickLabel',[], 'YTickLabel', [])
		set(h2,'XLim',get(h1,'XLim'),'Layer','top')
		legend(gca,state.imageProc.spineData.allDataList{:,1});
	end
end
set(gcf, 'Color' , [1 1 1]);
copyToClip(gcf);