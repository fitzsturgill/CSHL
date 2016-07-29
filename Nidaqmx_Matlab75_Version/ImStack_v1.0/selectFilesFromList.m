function names = selectFilesFromList(path, type)
global gh state

if nargin < 1
	path = pwd;
	type = '.tif'
end

if nargin == 1
	try
		filetype = get(gh.autotransformGUI.fileType, 'String');	
		value = get(gh.autotransformGUI.fileType, 'Value');
		filetype = filetype{value};
	catch
		filetype = type;
	end
end
if nargin == 2 
	filetype = type;
end

d = dir(fullfile(path, ['/*' filetype]));
if length(d) == 0
	str = 'No Files Found';
else
	str = {d.name};
end
str = sortrows({d.name}');
if isempty(str)
	names = {};
	disp('No Files with that Extension Selected. Please choose a Path with Image files');
	return
end

[s,v] = listdlg('PromptString','Select a file:', 'OKString', 'OK',...
	'SelectionMode','multiple',...
	'ListString', str, 'Name', 'Select a File');
names = str(s);
