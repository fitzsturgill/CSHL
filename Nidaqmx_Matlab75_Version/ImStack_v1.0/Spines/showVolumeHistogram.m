function showVolumeHistogram(type)
global gh state

if nargin < 1
	figure('name', ['Spine Volume Histogram for ' state.imageProc.spineData.spineDataName], 'NumberTitle', 'off');
	n = bar(state.imageProc.spineData.volumeHist(:,1), state.imageProc.spineData.volumeHist(:,2));
	set(get(gca,'YLabel'),'String','Number of Spines');
	set(get(gca,'XLabel'),'String','Spine Volume (fL)');
	h1 = gca;
	h2 = axes('Position',get(h1,'Position'));
	plot(state.imageProc.spineData.volumeHist(:,1),state.imageProc.spineData.volumeHist(:,2),'LineWidth',2)
	set(h2,'YAxisLocation','left','Color','none','XTickLabel',[],'YTickLabel', [])
	set(h2,'XLim',get(h1,'XLim'),'Layer','top')
	title( ['Spine Volume Histogram for ' state.imageProc.spineData.spineDataName]);

	
else
	if strcmp(type, 'all')
		figure('name', 'Spine Volume Histogram for All Files', 'NumberTitle', 'off');
		for i = 1:size(state.imageProc.spineData.allDataList,1)
			data(:,i) = state.imageProc.spineData.allDataList{i,4}(:,2);
		end
		n = bar(state.imageProc.spineData.volumeHist(:,1),data, 'group');
		set(get(gca,'YLabel'),'String','Number of Spines');
		set(get(gca,'XLabel'),'String','Spine Volume (fL)');
		title( ['Spine Volume Histogram for All Files']);
		h1 = gca;
		h2 = axes('Position',get(h1,'Position'));
		plot(state.imageProc.spineData.volumeHist(:,1),data,'LineWidth',2)
		set(h2,'YAxisLocation','left','Color','none','XTickLabel',[], 'YTickLabel', [])
		set(h2,'XLim',get(h1,'XLim'),'Layer','top')
		legend(gca,state.imageProc.spineData.allDataList{:,1});
	end
end
set(gcf, 'Color' , [1 1 1]);
copyToClip(gcf);