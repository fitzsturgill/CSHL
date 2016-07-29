function showOverlayImage
global gh state RGB

redchannel = 0;
bluechannel=0;
greenchannel=0;

for colorcounter = 1:3
	redchannelcounter = eval(['state.imageProc.overlay.red' num2str(colorcounter)]);
	if redchannelcounter
		redchannel = colorcounter;
		eval(['redFrame = state.imageProc.overlay.frameCounter' num2str(colorcounter) ';']);
	end
	
end

for colorcounter = 1:3
	greenchannelcounter = eval(['state.imageProc.overlay.green' num2str(colorcounter)]);
	if greenchannelcounter
		greenchannel = colorcounter;
		eval(['greenFrame = state.imageProc.overlay.frameCounter' num2str(colorcounter) ';']);
	end
end

for colorcounter = 1:3
	bluechannelcounter = eval(['state.imageProc.overlay.blue' num2str(colorcounter)]);
	if bluechannelcounter
		bluechannel = colorcounter;
		eval(['blueFrame = state.imageProc.overlay.frameCounter' num2str(colorcounter) ';']);
	end
end

if redchannel == 0
	redchannel = [];
end

if greenchannel == 0
	greenchannel = [];
end

if bluechannel == 0
	bluechannel = [];
end

if ~isempty(redchannel)
	eval(['valueRed = get(gh.overlayGUI.fileName' num2str(redchannel) ' , ''Value'');']);
	eval(['redImage = state.imageProc.overlay.weight' num2str(redchannel) '*double(state.imageProc.cell.currentImage{valueRed}(:,:,redFrame));']);
	rows = size(redImage,1);
	columns = size(redImage,2);
else
	redImage = [];
end

if ~isempty(greenchannel)
	eval(['valueGreen = get(gh.overlayGUI.fileName' num2str(greenchannel) ' , ''Value'');']);
	eval(['greenImage = state.imageProc.overlay.weight' num2str(greenchannel) '*double(state.imageProc.cell.currentImage{valueGreen}(:,:,greenFrame));']);
	rows = size(greenImage,1);
	columns = size(greenImage,2);
else
	greenImage = [];
end

if ~isempty(bluechannel)
	eval(['valueBlue = get(gh.overlayGUI.fileName' num2str(bluechannel) ' , ''Value'');']);
	eval(['blueImage = state.imageProc.overlay.weight' num2str(bluechannel) '*double(state.imageProc.cell.currentImage{valueBlue}(:,:,blueFrame));']);
	rows = size(blueImage,1);
	columns = size(blueImage,2);
else
	blueImage = [];
end


	if isempty(redImage)
		redImage = zeros(rows,columns);
	end

	if isempty(greenImage)
		greenImage = zeros(rows,columns);
	end
	
	if isempty(blueImage)
		blueImage = zeros(rows,columns);
	end
	
	RGB = double(cat(3,redImage,greenImage, blueImage));
	maxpixel = max(max(max(RGB)));
	RGB = 1/maxpixel*(RGB);
	RGB = (state.imageProc.overlay.overweight)*RGB;
	[a] = find(RGB > 1);
	RGB(a) = 1;
	r = size(RGB,1);
	c = size(RGB,2);
	state.imageProc.overlay.figure=figure('NumberTitle', 'off', 'Name', 'Overlayed Image', 'DoubleBuffer', 'on', 'pos', [300 300 c r]);
	rgbimagehandle = image('CData', RGB, 'CDataMapping', 'scaled');
	
	axis image
	set(gca, 'Clim', [0 1], 'YDir', 'Reverse','XTickLabelMode', 'manual', 'XTickLabel', [] , 'YTickLabel', [],'YTickLabelMode', 'manual', 'Position', [0 0 1 1]);
	




		
		
	