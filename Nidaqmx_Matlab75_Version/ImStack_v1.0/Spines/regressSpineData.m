function regressSpineData(type)
global gh state 

if nargin < 1
	type = 'denlen';
end

if isempty(state.imageProc.spineData.openPath)
	curDirSpineData;
else
	cd(state.imageProc.spineData.openPath);
end

filenames = selectFilesFromList(state.imageProc.spineData.openPath, '.txt');

if isempty(filenames)
	return
end

h = waitbar(0,'Starting Regression of Spine Data Files....', 'Name', 'Regression of Spine Statistics', 'Pointer', 'watch', ...
	'Position', [483.750000000000   305.250000000000   270.000000000000   056.250000000000]);
total = length(filenames);
% Make large arrays, one for each parameter in each file.
for i = 1:total
	waitbar((i-1)/total,h, filenames{i});
	[annovaName, zerosRemovedAnova,  anDensity2d, anDensity3d, anLen, anVol]  = anRegSpineDensityVar([state.imageProc.spineData.openPath filenames{i}]);
	name = filenames{i}(1:end-4);
	switch type
	case 'denlen'
		figure('NumberTitle', 'off');
		regS = [ones(size(anLen,2),1) anLen'];
		size(anDensity2d')
		size(regS)
		[a,b,c,d,r]  = regress(anDensity2d', regS);
		maxLen = max(anLen);
		maxden = max(anDensity2d);
		x = linspace(0,maxLen, 100);
		y = a(1)*x+ a(2);
		plot(anLen,anDensity2d,'+',x,y);
		set(get(gca, 'YLabel'),'String', '2D Density (per um)');
		set(get(gca, 'XLabel'), 'String', 'Mean Spine Length (um)');
		set(gca, 'XLim', [0 (1.1*maxLen)], 'YLim', [ 0 1.1*maxden]);
		title(['R^2 = ' num2str(r(1)) ' for Spine 2D Density vs. Length for ' name]);
		text(.8*maxLen,.8*maxden,['Y = ' num2str(a(1)) ' X + ' num2str(a(2))]);

	case 'lenvol'
		figure('NumberTitle', 'off');
		regS = [ones(size(anLen,2),1) anLen'];
		[a,b,c,d,r] = regress(anVol',regS);
		maxLen = max(anLen);
		maxVol = max(anVol);
		x = linspace(0,maxLen, 100);
		y = a(1)*x+a(2);
		plot(anLen, anVol, '+', x,y);
		set(get(gca, 'YLabel'), 'String', 'Mean Spine Volume (fL)');
		set(get(gca, 'XLabel'), 'String', 'Mean Spine Length (um)');
		set(gca, 'XLim', [0 (1.1*maxLen)], 'YLim', [ 0 1.1*maxVol]);
		title(['R^2 = ' num2str(r(1)) ' for Spine Volume vs. Length for ' name]);
		text(.8*maxLen,.8*maxVol,['Y = ' num2str(a(1)) ' X + ' num2str(a(2))]);

	case 'denvol'
		figure('NumberTitle', 'off');
		regS = [ones(size(anVol,2),1) anVol'];
		[a,b,c,d,r]  = regress(anDensity2d', regS);
		maxVol = max(anVol);
		maxden = max(anDensity2d);
		x = linspace(0,maxVol, 100);
		y = a(1)*x+ a(2);
		plot(anVol,anDensity2d,'+',x,y);
		set(get(gca, 'YLabel'),'String', '2D Density (per um)');
		set(get(gca, 'XLabel'), 'String', 'Mean Spine Volume (fL)');
		set(gca, 'XLim', [0 (1.1*maxVol)], 'YLim', [ 0 1.1*maxden]);
		title(['R^2 = ' num2str(r(1)) ' for Spine 2D Density vs. Volume for ' name]);
		text(.8*maxVol,.8*maxden,['Y = ' num2str(a(1)) ' X + ' num2str(a(2))]);

	case 'all'
		figure('NumberTitle', 'off');
		regS = [ones(size(anLen,2),1) anLen'];
		[a,b,c,d,r]  = regress(anDensity2d', regS);
		maxLen = max(anLen);
		maxden = max(anDensity2d);
		x = linspace(0,maxLen, 100);
		y = a(1)*x+ a(2);
		plot(anLen,anDensity2d,'+',x,y);
		set(get(gca, 'YLabel'),'String', '2D Density (per um)');
		set(get(gca, 'XLabel'), 'String', 'Mean Spine Length (um)');
		set(gca, 'XLim', [0 (1.1*maxLen)], 'YLim', [ 0 1.1*maxden]);
		title(['R^2 = ' num2str(r(1)) ' for Spine 2D Density vs. Length for ' name]);
		text(.8*maxLen,.8*maxden,['Y = ' num2str(a(1)) ' X + ' num2str(a(2))]);

		figure('NumberTitle', 'off');
		regS = [ones(size(anLen,2),1) anLen'];
		[a,b,c,d,r] = regress(anVol',regS);
		maxLen = max(anLen);
		maxVol = max(anVol);
		x = linspace(0,maxLen, 100);
		y = a(1)*x+a(2);
		plot(anLen, anVol, '+', x,y);
		set(get(gca, 'YLabel'), 'String', 'Mean Spine Volume (fL)');
		set(get(gca, 'XLabel'), 'String', 'Mean Spine Length (um)');
		set(gca, 'XLim', [0 (1.1*maxLen)], 'YLim', [ 0 1.1*maxVol]);
		title(['R^2 = ' num2str(r(1)) ' for Spine Volume vs. Length for ' name]);
		text(.8*maxLen,.8*maxVol,['Y = ' num2str(a(1)) ' X + ' num2str(a(2))]);

		figure('NumberTitle', 'off');
		regS = [ones(size(anVol,2),1) anVol'];
		[a,b,c,d,r]  = regress(anDensity2d', regS);
		maxVol = max(anVol);
		maxden = max(anDensity2d);
		x = linspace(0,maxVol, 100);
		y = a(1)*x+ a(2);
		plot(anVol,anDensity2d,'+',x,y);
		set(get(gca, 'YLabel'),'String', '2D Density (per um)');
		set(get(gca, 'XLabel'), 'String', 'Mean Spine Volume (fL)');
		set(gca, 'XLim', [0 (1.1*maxVol)], 'YLim', [ 0 1.1*maxden]);
		title(['R^2 = ' num2str(r(1)) ' for Spine 2D Density vs. Volume for ' name]);
		text(.8*maxVol,.8*maxden,['Y = ' num2str(a(1)) ' X + ' num2str(a(2))]);
	end
	waitbar(i/total,h);
end
close(h);





