function [totalframes, pixelsPerLine, linesPerFrame, numberOfZSlicesPerChannel] = extractImageInfoFromParsedValues(image, ...
		numberOfChannels, average, framesPerSlice)
global gh state

% this function will determine the format of the image from the parsed values.

pixelsPerLine = size(image,2);
linesPerFrame = size(image,1);
totalframes = size(image,3);

numberOfZSlicesPerChannel = totalframes/(framesPerSlice*numberOfChannels);
	