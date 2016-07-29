function out = doMultiAnova(type, compare, display)
global gh state

% parses a text file from 3DMA and does a comparison of variations.  Output is the p value sfrom the Anova F test
% It will also compare graphically the sorces of variations.

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


[fname, pname] = uigetfile('*.txt', 'Choose Data to ANOVA...');
if fname < 1
	return
else
	filename = [pname fname];
	state.imageProc.spineData.openPath = pname;
end	
set(gh.spineDataGUI.figure1, 'pointer', 'watch');
cd(pname);

[annovaName, zerosRemovedAnova,  anDensity2d, anDensity3d, anLen, anVol]  = anSpineDensVar(filename);

switch type
case 'density2d'
	[p2d,tbl2d,stat2d] = anova1(anDensity2d, annovaName, display);
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
	[p3d,tbl3d,stat3d]= anova1(anDensity3d, annovaName,  display);
	if strcmp(compare, 'yes')
		figure('pos', [232   130   661   548]);
		[c3d,m3d] = multcompare(stat3d);
	end
	pairs = showSignificantPairs(real(c3d),stat3d, '3D Density');
	title('3D Density Comparison');
	setLineColors(gca,'black');
	copyToClip(gcf);
	out = p3d
case 'length'
	[pLen,tblLen,statLen]= anova1(anLen, zerosRemovedAnova, display);
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
	[pVol,tblVol,statVol] = anova1(anVol, zerosRemovedAnova,  display);
	if strcmp(compare, 'yes')
		figure('pos', [232   130   661   548]);
		[cVol,mVol] = multcompare(statVol);
	end
	pairs = showSignificantPairs(real(cVol),statVol, 'Volume');
	title('Volume Comparison');
	setLineColors(gca,'black');
	copyToClip(gcf);
	out = pVol;
case 'all'
	[p2d,tbl2d,stat2d] = anova1(anDensity2d, annovaName,display);
	[p3d,tbl3d,stat3d]= anova1(anDensity3d, annovaName, display);
	[pLen,tblLen,statLen]= anova1(anLen, zerosRemovedAnova, display);
	[pVol,tblVol,statVol] = anova1(anVol, zerosRemovedAnova,  display);
	if strcmp(compare, 'yes')
		figure('pos', [232   130   661   548]);
		[c2d,m2d] = multcompare(stat2d);
		title('2D Density Comparison');
		setLineColors(gca,'black');
		figure('pos', [232   130   661   548]);
		[c3d,m3d] = multcompare(stat3d);
		title('3D Density Comparison');
		setLineColors(gca,'black');
		figure('pos', [232   130   661   548]);
		[cLen,mLen] = multcompare(statLen);
		title('Length Comparison');
		setLineColors(gca,'black');
		figure('pos', [232   130   661   548]);
		[cVol,mVol] = multcompare(statVol);
		title('Volume Comparison');
		setLineColors(gca,'black');
		pairs2d = showSignificantPairs(real(c2d),stat2d, '2D Density');
		pairs3d = showSignificantPairs(real(c3d),stat3d, '3D Density');
		pairsLen = showSignificantPairs(real(cLen),statLen, 'Length');
		pairsVol = showSignificantPairs(real(cVol),statVol, 'Volume');
	end
	out = [p2d p3d pLen pVol];
end

set(gh.spineDataGUI.figure1, 'pointer', 'arrow');