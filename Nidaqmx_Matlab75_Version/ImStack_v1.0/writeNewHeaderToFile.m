function writeNewHeaderToFile(image, header)
global gh state

data = openTif(image);
[path,name,ext] = fileparts(image);

if strcmp(path, '') % no path,. use current
	arrayToTiff(data, name, header);
else
	arrayToTiff(data, [path '\' name], header);
end
