function out = doOverallKrusWal(type, compare, display)
global gh state 

if nargin < 1
	type = 'density2d';
	compare = 'yes';
	display = 'off';
elseif nargin < 2
	compare = 'yes';
	display = 'off';
elseif nargin < 3
	display = 'off';
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

h = waitbar(0,'Starting To Kruskal-Wallis Spine Data Files....', 'Name', 'Kruskal-Wallis of Constructs', 'Pointer', 'watch', ...
	'Position', [483.750000000000   305.250000000000   270.000000000000   056.250000000000]);
total = length(filenames);

total2dDensity = [];
total3dDensity = [];
totalLen = [];
totalVol = [];
totalAnnovaZerName = {};
totalAnnovaName = {};

% Make large arrays, one for each parameter in each file.
for i = 1:total
	waitbar((i-1)/total,h, filenames{i});
	[annovaName, zerosRemovedAnova,  anDensity2d, anDensity3d, anLen, anVol]  = anSpineDensVar([state.imageProc.spineData.openPath filenames{i}]);
	name = filenames{i}(1:end-4);
	for j = 1:length(annovaName)
		tempAnnovaName{j} = name;
	end
	for j = 1:length(zerosRemovedAnova)
		tempAnnovaZerName{j} = name;
	end
	totalAnnovaName = [totalAnnovaName tempAnnovaName];
	totalAnnovaZerName = [totalAnnovaZerName tempAnnovaZerName];
	
	total2dDensity = [total2dDensity anDensity2d];
	total3dDensity = [total3dDensity anDensity3d];
	totalLen = [totalLen anLen];
	totalVol = [totalVol anVol];
	waitbar(i/total,h);
end
close(h);

switch type
case 'density2d'
	[p2d,tbl2d,stat2d] = kruskalwallis(total2dDensity, totalAnnovaName, display);
	if strcmp(compare, 'yes')
		figure('pos', [232   130   661   548]);
		[c2d,m2d] = multcompare(stat2d);
	end
	pairs = showSignificantPairs(real(c2d),stat2d, '2D Density');
	title('2D Density Comparison');
	setLineColors(gca,'black');
	copyToClip(gcf);
	out = p2d;
case 'density3d'
	[p3d,tbl3d,stat3d]= kruskalwallis(total3dDensity, totalAnnovaName,  display);
	if strcmp(compare, 'yes')
		figure('pos', [232   130   661   548]);
		[c3d,m3d] = multcompare(stat3d);
	end
	pairs = showSignificantPairs(real(c3d),stat3d, '3D Density');
	title('3D Density Comparison');
	setLineColors(gca,'black');
	copyToClip(gcf);
	out = p3d;
case 'length'
	[pLen,tblLen,statLen]= kruskalwallis(totalLen, totalAnnovaZerName, display);
	if strcmp(compare, 'yes')
		figure('pos', [232   130   661   548]);
		[cLen,mLen] = multcompare(statLen);
	end
	pairs = showSignificantPairs(real(cLen),statLen, 'Length');
	title('Length Comparison');
	setLineColors(gca,'black');
	copyToClip(gcf);
	out = pLen;
case 'volume'
	[pVol,tblVol,statVol] = kruskalwallis(totalVol, totalAnnovaZerName,  display);
	if strcmp(compare, 'yes')
		figure('pos', [232   130   661   548]);
		[cVol,mVol] = multcompare(statVol);
	end
	pairs = showSignificantPairs(real(cVol), statVol, 'Volume');
	title('Volume Comparison');
	setLineColors(gca,'black');
	copyToClip(gcf);
	out = pVol;
case 'all'
	[p2d,tbl2d,stat2d] = kruskalwallis(total2dDensity, totalAnnovaName, display);
	if strcmp(compare, 'yes')
		figure('pos', [232   130   661   548]);
		[c2d,m2d] = multcompare(stat2d);
	end
	pairs = showSignificantPairs(real(c2d),stat2d, '2D Density');
	title('2D Density Comparison');
	setLineColors(gca,'black');
	copyToClip(gcf);
	[p3d,tbl3d,stat3d]= kruskalwallis(total3dDensity, totalAnnovaName,  display);
	if strcmp(compare, 'yes')
		figure('pos', [232   130   661   548]);
		[c3d,m3d] = multcompare(stat3d);
	end
	pairs = showSignificantPairs(real(c3d),stat3d, '3D Density');
	title('3D Density Comparison');
	setLineColors(gca,'black');
	copyToClip(gcf);
	
	[pLen,tblLen,statLen]= kruskalwallis(totalLen, totalAnnovaZerName, display);
	if strcmp(compare, 'yes')
		figure('pos', [232   130   661   548]);
		[cLen,mLen] = multcompare(statLen);
	end
	pairs = showSignificantPairs(real(cLen),statLen, 'Length');
	title('Length Comparison');
	setLineColors(gca,'black');
	copyToClip(gcf);
	out = pLen;
	
	[pVol,tblVol,statVol] = kruskalwallis(totalVol, totalAnnovaZerName,  display);
	if strcmp(compare, 'yes')
		figure('pos', [232   130   661   548]);
		[cVol,mVol] = multcompare(statVol);
	end
	pairs = showSignificantPairs(real(cVol),statVol, 'Volume');
	title('Volume Comparison');
	setLineColors(gca,'black');
	copyToClip(gcf);
	out = [p2d p3d pLen pVol];
	
end



	

