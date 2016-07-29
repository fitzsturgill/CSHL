function newimage = extractData(image, numberofchannels)
global gh state

%function breaks the image up into a cell array where the cell 
% index is the channel number

if numberofchannels == 1
	newimage{1} = image;
	return
end
totalframes = size(image,3);
if state.imageProc.internal.interleaved
	for channelcounter = 1:numberofchannels
		newimagecounter = 1;
		for framecounter = 0:(totalframes/numberofchannels-1)
			newimage{channelcounter}(:,:,newimagecounter) = ...
				image(:,:,(channelcounter + framecounter*numberofchannels));
			newimagecounter = newimagecounter+1;
		end
	end
else
	totalFrames = size(image,3);
	number = 0;
	for i = 1:numberofchannels
		newimage{i} = image(:,:,(1 + (i-1)*totalFrames/numberofchannels):(i*totalFrames/numberofchannels));
	end
end
