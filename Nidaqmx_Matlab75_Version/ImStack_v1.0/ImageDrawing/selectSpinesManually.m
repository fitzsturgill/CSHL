function selectSpinesManually(axis)
global gh state bit
% bit is 0 if it is a normal point and 1 if it si an extension.
bit = [];
peaks=[];

XLim = get(axis, 'XLim');
YLim = get(axis, 'YLim');

[X,Y] = getptsSpine(axis);
totalPoints = length(X);
if totalPoints <=1
	return
end
pointPositions=round([X Y]);
pointPositions = [pointPositions bit'];
numberOfSpines = floor((size(pointPositions,1) - sum(bit))/2);
oldNSpines=state.imageProc.spine.numberOfSpines;

state.imageProc.spine.numberOfSpines = state.imageProc.spine.numberOfSpines + numberOfSpines;
updateGUIByGlobal('state.imageProc.spine.numberOfSpines');

val=get(gh.spineGUI.spineColorMap, 'Value');
str=get(gh.spineGUI.spineColorMap, 'String');

str = str{val};
switch str
case 'Jet'
	map = jet(numberOfSpines);
case 'Green'
	map = repmat([0 1 0],numberOfSpines,1);
case 'Blue'
	map = repmat([0 0 1],numberOfSpines,1);
case 'Red '
	map = repmat([1 0 0],numberOfSpines,1);
case 'Gray'
	map = gray(numberOfSpines);
otherwise
	map = jet(numberOfSpines);
end



spacer = 0;
skip = 0;
for lineCounter = 1:numberOfSpines
	
	Xpoints = pointPositions((lineCounter + spacer):(lineCounter + spacer+1),1);
	Ypoints = pointPositions((lineCounter + spacer):(lineCounter + spacer+1),2);
	
	for next = 1:3
		if (lineCounter + spacer+ 1) < totalPoints
			try
				if pointPositions((lineCounter + spacer + 1 + next),3) 
					Xpoints = [Xpoints; pointPositions((lineCounter + spacer + 1 + next ),1)];
					Ypoints = [Ypoints; pointPositions((lineCounter + spacer + 1 + next ),2)];
					skip = skip+1;
				else
					break
				end
			end
		end
	end
	% draw connecting Lines
	spacer = spacer+skip+1;
	skip=0;
	state.imageProc.spine.linehandles(lineCounter+oldNSpines) = line(Xpoints, Ypoints, 'tag', ['Spine ' num2str(lineCounter+oldNSpines)], 'color', map(lineCounter,:), ...
		'linewidth', state.imageProc.spine.spineLine,'Parent',gh.spineGUI.mainAxes);
	state.imageProc.spine.texthandles(lineCounter+oldNSpines)= text(Xpoints(end)+4, Ypoints(end)+4, [num2str(lineCounter+oldNSpines)], 'color', map(lineCounter,:), 'fontsize', ...
		state.imageProc.spine.spineText,'Parent',gh.spineGUI.mainAxes);
	if state.imageProc.spine.calcSpineVols
		% Calc Vols?
		[int]=improfile(state.imageProc.spine.maniImgData,Xpoints, Ypoints);
		maxBinary=imextendedmax(int,10);	%Finds maxima
		NewPeaks=double(maxBinary).*double(int)
		NewPeaks=NewPeaks(imregionalmax(NewPeaks))%	Get main peaks out...
		maxVal=1;
		maxVal2=1;
		if length(NewPeaks) > 2
% 			NewPeaks=imhmax(NewPeaks,200);
% 			NewPeaks=NewPeaks(imregionalmax(NewPeaks));
			NewPeaks = unique(NewPeaks);
			maxVal=NewPeaks(1);
			if length(NewPeaks) >= 2
				maxVal2=NewPeaks(2);
			else
				maxVal2=maxVal;
			end
		elseif length(NewPeaks) == 2
			maxVal=NewPeaks(1);
			maxVal2=NewPeaks(2);
		elseif length(NewPeaks) == 1
			maxVal=NewPeaks(1);
			maxVal2=maxVal;
		end
		Ratio=maxVal2/maxVal;
		if lineCounter==1
			f=figure('color','w','DoubleBuffer','on','NumberTitle','off','Name','Pixel Intensitiues along spines');
			subplot(numberOfSpines,1,lineCounter); 
			plot(int); 
			title(['Intensity accross Spine number ' ...
					num2str(lineCounter) ' Ratio of Peaks (second/first) = ' num2str(Ratio) '.']);
		else
			set(f,'visible','on');
			subplot(numberOfSpines,1,lineCounter); 
			plot(int); 
			title(['Intensity accross Spine number ' ...
					num2str(lineCounter) ' Ratio of Peaks (second/first) = ' num2str(Ratio) '.']);
		end
		peaks=[peaks maxVal maxVal2];
	end
end

set(axis, 'XLim', XLim, 'YLim', YLim);
state.imageProc.spine.spineLengths = computeLengthOfLine(state.imageProc.spine.linehandles);
columns = size(state.imageProc.spine.spineLengths,2);
state.imageProc.spine.spineDensity = columns/state.imageProc.spine.dendriteLength;
updateGUIByGlobal('state.imageProc.spine.spineDensity');

saveSpineDataToExcel;