function header = readImageHeaderTif(imagename)
global state gh

try
	a = imfinfo(imagename);
	sizeOfHeader = length(a);
	if sizeOfHeader > 1
		header = a(1).ImageDescription;
	else
		header = a.ImageDescription;
	end
	
catch
	header = 'No Image Description Found.';
end
