function compareSpineStats(stat)

if nargin < 1
	stat = 'density';
end

global gh state
flag = 1;
for i = 1:size(state.imageProc.spineData.allDataList,1)	
	names{i} = state.imageProc.spineData.allDataList{i,1};
end
for i = 1:size(state.imageProc.spineData.allDataList,1)	
	if state.imageProc.spineData.allDataList{i,2}(end) > 1 & state.imageProc.spineData.analyzeByNeuron == 1
		label{i} = ['(N = ' num2str(state.imageProc.spineData.allDataList{i,2}(end)) ')'];
	else
		label = '';
	end
end

switch stat
case 'Length'
	stat = 'Length (um)';
	for i = 1:size(state.imageProc.spineData.allDataList,1)
		data(i) = state.imageProc.spineData.allDataList{i,2}(7);
		err(i) = state.imageProc.spineData.allDataList{i,2}(8);
	end
	data = data';
	err = err';
case 'Density'
	stat = '2D Density (per um)';
	for i = 1:size(state.imageProc.spineData.allDataList,1)
		data(i) = state.imageProc.spineData.allDataList{i,2}(1);
		err(i) = state.imageProc.spineData.allDataList{i,2}(2);
	end
	data = data';
	err = err';
case 'Density3'
	stat = '3D Density (per um)';
	for i = 1:size(state.imageProc.spineData.allDataList,1)
		data(i) = state.imageProc.spineData.allDataList{i,2}(13);
		err(i) = state.imageProc.spineData.allDataList{i,2}(14);
	end
	data = data';
	err = err';

case 'Volume'
	stat = 'Volume (fL)';
	for i = 1:size(state.imageProc.spineData.allDataList,1)
		data(i) = state.imageProc.spineData.allDataList{i,2}(11);
		err(i) = state.imageProc.spineData.allDataList{i,2}(12);
	end
	data = data';
	err = err';
case 'all'
	figure('name', 'Comparison of Spine Density, Length, and Volume', 'NumberTitle', 'off');
	for i = 1:size(state.imageProc.spineData.allDataList,1)
		data(i) = state.imageProc.spineData.allDataList{i,2}(7);
		err(i) = state.imageProc.spineData.allDataList{i,2}(8);
	end
	lendata = data';
	lenerr = err';
	
	for i = 1:size(state.imageProc.spineData.allDataList,1)
		data(i) = state.imageProc.spineData.allDataList{i,2}(11);
		err(i) = state.imageProc.spineData.allDataList{i,2}(12);
	end
	voldata = data';
	volerr = err';
	
	for i = 1:size(state.imageProc.spineData.allDataList,1)
		data(i) = state.imageProc.spineData.allDataList{i,2}(1);
		err(i) = state.imageProc.spineData.allDataList{i,2}(2);
	end
	dendata = data';
	denerr = err';
	
	for i = 1:size(state.imageProc.spineData.allDataList,1)
		data(i) = state.imageProc.spineData.allDataList{i,2}(13);
		err(i) = state.imageProc.spineData.allDataList{i,2}(14);
	end
	dendata3 = data';
	denerr3 = err';
	
	subplot(1,4,1);
	plotBarAndErr(dendata,denerr);
	set(get(gca,'YLabel'),'String','2D Density (per um)');
	set(get(gca,'XLabel'),'String','Name of Construct');
	set(gca,'XTickLabel', names);
	set(gca, 'YLim', [0 1.15*(max(dendata)+max(denerr))]);
	title( ['2D Density']);
	
	subplot(1,4,2);
	plotBarAndErr(dendata3,denerr3);
	set(get(gca,'YLabel'),'String','3D Density (per um)');
	set(get(gca,'XLabel'),'String','Name of Construct');
	set(gca,'XTickLabel', names);
	set(gca, 'YLim', [0 1.15*(max(dendata)+max(denerr))]);	
	title( ['3D Density']);

	
	subplot(1,4,3);
	plotBarAndErr(lendata,lenerr);
	set(get(gca,'YLabel'),'String','Length (um)');
	set(get(gca,'XLabel'),'String','Name of Construct');
	set(gca,'XTickLabel', names);
	set(gca, 'YLim', [0 1.15*(max(lendata)+max(lenerr))]);
	title( ['Mean Length']);

	
	subplot(1,4,4);
	plotBarAndErr(voldata,volerr);
	set(get(gca,'YLabel'),'String','Volume (fL)');
	set(get(gca,'XLabel'),'String','Name of Construct');
	set(gca,'XTickLabel', names);
	set(gca, 'YLim', [0 1.15*(max(voldata)+max(volerr))]);
	title( ['Mean Volume']);


	flag = 0;
end

if flag
	figure('name', ['Comparison of Spine ' stat ], 'NumberTitle', 'off');
	center = plotBarAndErr(data,err);
	title( ['Comparison of Spine ' stat ]);
	set(get(gca,'YLabel'),'String',stat);
	set(get(gca,'XLabel'),'String','Name of Construct');
	set(gca,'XTickLabel', names);
	set(gca, 'YLim', [0 1.15*(max(data)+max(err))]);
	height = 1.05*(data+err);
	text(center, height, label, 'HorizontalAlignment','center');
end

set(gcf, 'Color' , [1 1 1]);
copyToClip(gcf);

