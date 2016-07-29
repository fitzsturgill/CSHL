function saveAndLoadExtractedData %(image, numberofchannels, filename, comments)
global gh state

filename = get(gh.imageProcessingGUI.fileName, 'String');
value = get(gh.imageProcessingGUI.fileName, 'Value');
if ischar(filename)
	[path,name,ext,ver] = fileparts(filename);
	cd(path);
else
	[path,name,ext,ver] = fileparts(filename{value});
	cd(path);
end

image = state.imageProc.cell.currentImage{value};
numberofchannels = state.imageProc.parsing.cell.numberOfChannels{value};
comments = '';

% this function will extract the data from an image and write it as a tif file.

		
% creates cell array with extracted data.

	newimage = extractData(image, numberofchannels);
	
	for channelcounter = 1:numberofchannels
		[fname, pname] = uiputfile('*.tif', ['Save Extracted Channel Number ' ...
			num2str(channelcounter) ' Image as...' ]);
		if pname > 0
			
			newfilename = [pname fname];
			arrayToTiff(newimage{channelcounter}, newfilename, comments);
			loadImageFromName([newfilename '.tif'], pname);
		else
			return
		end
	end
	
	
	