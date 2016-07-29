function closeAllImages
global gh state

string = get(gh.imageProcessingGUI.fileName, 'String');
if ischar(string)
	filesOpen = 1;
else
	filesOpen = length(string);
end
for i = 1:filesOpen
	set(gh.imageProcessingGUI.fileName, 'Value', 1);
	closeCurrentImage;
end
