function figureButtonDownFcn
global state gh

% This is the button down function fo rthe figure...iif pressed, it updates the GUI to make
% The info current to the image

name=get(gcf, 'Name');
value = getMenuIndex2(gh.imageProcessingGUI.fileName, name);
currentvalue = get(gh.imageProcessingGUI.fileName, 'String');


if strcmp(value, currentvalue)
	return
else
	switchFileName(value);
	set(gh.imageProcessingGUI.fileName, 'Value', value);
end

	