function extractChannelOfDataFromImage(image,numberOfChannels,average,framesPerSlice)
global state gh

if channelNumber == 1
	return
else
	[totalframes, pixelsPerLine, linesPerFrame, numberOfZSlicesPerChannel] = extractImageInfoFromParsedValues(image, ...
		numberOfChannels, average, framesPerSlice);
	