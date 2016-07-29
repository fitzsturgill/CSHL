function [text, densityNumber, densityError, overdensityNumber, overSpineLength, overSpineLenError, ...
		meanSpineLength, meanSpineLenError, overSpineVolume, overSpineVolError, ...
		meanSpineVolume, meanSpineVolError, lengthHist, volHist] = readSpineData(filename)
global gh state 

if state.imageProc.spineData.analyzeByFile==1
	
	text = textread(filename, '%s', 'headerlines', 1, 'delimiter', '\n');
	
	counter = 0;
	for i = 1:length(text)
		index = findstr('Mean spine density', text{i});
		if ~isempty(index)
			counter = counter + 1;
			density{counter} = text{i};
			tokens = tokenize(density{counter});
			for j = 1:length(tokens)
				place = findstr('density', tokens{j});
				if ~isempty(place)
					densityNumber(counter) = str2num(tokens{j+1});
					densityError(counter) = str2num(tokens{j+3});
				end
			end
		end	
	end
	
	
	counter = 0;
	for i = 1:length(text)
		index = findstr('Overall spine density (Total #spine/Total length of dendrites)' , text{i});
		if ~isempty(index)
			counter = counter + 1;
			density{counter} = text{i};
			tokens = tokenize(density{counter});
			for j = 1:length(tokens)
				place = findstr('dendrites)', tokens{j});
				if ~isempty(place)
					try
						overdensityNumber(counter) = str2num(tokens{j+1});
					end		
				end
			end
		end	
	end
	
	counter = 0;
	for i = 1:length(text)
		index = findstr('Sample mean of spine lengths ' , text{i});
		if ~isempty(index)
			counter = counter + 1;
			density{counter} = text{i};
			tokens = tokenize(density{counter});
			for j = 1:length(tokens)
				place = findstr('lengths', tokens{j});
				if ~isempty(place)
					meanSpineLength(counter) = str2num(tokens{j+1});
					meanSpineLenError(counter) = str2num(tokens{j+3});
				end
			end
		end	
	end
	
	counter = 0;
	for i = 1:length(text)
		index = findstr('Overall average spine lengths ' , text{i});
		if ~isempty(index)
			counter = counter + 1;
			density{counter} = text{i};
			tokens = tokenize(density{counter});
			for j = 1:length(tokens)
				place = findstr('lengths', tokens{j});
				if ~isempty(place)
					overSpineLength(counter) = str2num(tokens{j+1});
					overSpineLenError(counter) = str2num(tokens{j+3});
				end
			end
		end	
	end
	
	counter = 0;
	for i = 1:length(text)
		index = findstr('Sample mean of spine volumes ' , text{i});
		if ~isempty(index)
			counter = counter + 1;
			density{counter} = text{i};
			tokens = tokenize(density{counter});
			for j = 1:length(tokens)
				place = findstr('volumes', tokens{j});
				if ~isempty(place)
					meanSpineVolume(counter) = str2num(tokens{j+1});
					meanSpineVolError(counter) = str2num(tokens{j+3});
				end
			end
		end	
	end
	
	counter = 0;
	for i = 1:length(text)
		index = findstr('Overall average spine volumes ' , text{i});
		if ~isempty(index)
			counter = counter + 1;
			density{counter} = text{i};
			tokens = tokenize(density{counter});
			for j = 1:length(tokens)
				place = findstr('volumes', tokens{j});
				if ~isempty(place)
					overSpineVolume(counter) = str2num(tokens{j+1});
					overSpineVolError(counter) = str2num(tokens{j+3});
				end
			end
		end	
	end
	
	counter = 1;
	for i = 1:length(text)
		index = findstr('Spine length distribution' , text{i});
		if ~isempty(index)
			startOfLenHist = i;
			for j = i:length(text)
				try
					lengthHist{counter} = text{i+counter};
					counter = counter+ 1;
				end
			end
		end	
	end
	
	for i = 1:26
		tokens = tokenize(lengthHist{i});
		bins(i) = str2num(tokens{1});
		number(i) = str2num(tokens{2});
	end
	
	lengthHist = [bins' number'];

	% Where does the spine Volume list end
	for i = 1:length(text)
		index = findstr('Sample mean of spine volumes ', text{i});
		if ~isempty(index)
			endOfVol = i;
			break
		end
	end

	counter = 1;
	for i = 1:length(text)
		index = findstr('Dendrite	Mean volume' , text{i});
		if ~isempty(index)
			for j = i+1:endOfVol-1
				tk = tokenize(text{j});
				volData(counter)= str2num(tk{2});
				counter = counter+ 1;
			end
		end	
	end
	minVal = min(volData);
	maxVal = max(volData);
	x = [ 0 linspace(minVal,maxVal, 26)];
	volHistogram = hist(volData,x);
	volHist = [x' volHistogram'];
	
	for i = 1:length(text)
		index = findstr('Total no. of spines:'  , text{i});
		if ~isempty(index)
			tokens = tokenize(text{i});
			state.imageProc.spineData.numberOfSpines = str2num(tokens{5});
			updateGUIByGlobal('state.imageProc.spineData.numberOfSpines');
			break
			
		end	
	end
	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif state.imageProc.spineData.analyzeByNeuron == 1
	
	text = textread(filename, '%s', 'headerlines', 1, 'delimiter', '\n');
	
	counter = 0;
	for i = 1:length(text)
		index = findstr('Overall average spine volumes ' , text{i});
		if ~isempty(index)
			counter = counter + 1;
			density{counter} = text{i};
			tokens = tokenize(density{counter});
			for j = 1:length(tokens)
				place = findstr('volumes', tokens{j});
				if ~isempty(place)
					overSpineVolume(counter) = str2num(tokens{j+1});
					overSpineVolError(counter) = str2num(tokens{j+3});
				end
			end
		end	
	end
	
	counter = 0;
	for i = 1:length(text)
		index = findstr('Overall spine density (Total #spine/Total length of dendrites)' , text{i});
		if ~isempty(index)
			counter = counter + 1;
			density{counter} = text{i};
			tokens = tokenize(density{counter});
			for j = 1:length(tokens)
				place = findstr('dendrites)', tokens{j});
				if ~isempty(place)
					try
						overdensityNumber(counter) = str2num(tokens{j+1});
					end		
				end
			end
		end	
	end
	
	counter = 0;
	for i = 1:length(text)
		index = findstr('Overall average spine lengths ' , text{i});
		if ~isempty(index)
			counter = counter + 1;
			density{counter} = text{i};
			tokens = tokenize(density{counter});
			for j = 1:length(tokens)
				place = findstr('lengths', tokens{j});
				if ~isempty(place)
					overSpineLength(counter) = str2num(tokens{j+1});
					overSpineLenError(counter) = str2num(tokens{j+3});
				end
			end
		end	
	end
	
	counter = 1;
	for i = 1:length(text)
		index = findstr('Spine length distribution' , text{i});
		if ~isempty(index)
			startOfLenHist = i;
			for j = i:length(text)
				try
					lengthHist{counter} = text{i+counter};
					counter = counter+ 1;
				end
			end
		end	
	end

	for i = 1:26
		tokens = tokenize(lengthHist{i});
		bins(i) = str2num(tokens{1});
		number(i) = str2num(tokens{2});
	end
	
	lengthHist = [bins' number'];
	% Where does the spine Volume list end
	for i = 1:length(text)
		index = findstr('Sample mean of spine volumes ', text{i});
		if ~isempty(index)
			endOfVol = i;
			break
		end
	end

	counter = 1;
	for i = 1:length(text)
		index = findstr('Dendrite	Mean volume' , text{i});
		if ~isempty(index)
			for j = i+1:endOfVol-1
				tk = tokenize(text{j});
				volData(counter)= str2num(tk{2});
				counter = counter+ 1;
			end
		end	
	end
	minVal = min(volData);
	maxVal = max(volData);
	x = [ 0 linspace(minVal,maxVal, 26)];
	volHistogram = hist(volData,x);
	volHist = [x' volHistogram'];

	for i = 1:length(text)
		index = findstr('Total no. of spines:'  , text{i});
		if ~isempty(index)
			tokens = tokenize(text{i});
			state.imageProc.spineData.numberOfSpines = str2num(tokens{5});
			updateGUIByGlobal('state.imageProc.spineData.numberOfSpines');
			break
		end	
	end
	
% Break the data apart per neurons	
	counter = 0;
	for i = 1:length(text)
		index = findstr('Dendrite	#spine', text{i});
		if ~isempty(index)
			start = i +1;
			for k = i:length(text)
				endOfData = findstr('Dendrite	Spine Density(2D)', text{k});
				if ~isempty(endOfData)
					endO = k-2;
					break
				end
			end
			break
		end
	end
	
	totalFiles = endO-start+1;
	dencount =1;
	zercount = 1;
	name = 'test';
	for j = start:endO
		tokens = tokenize(text{j});
		comma = findstr(',', tokens{5})+1;
		if strcmp(name, tokens{5}(comma:end-12)) 
			
			if strcmp(tokens{2}, '0')
				zer(zercount) = j;
				zercount = zercount+1;	
			end
		else
			newDendrite(dencount) = j;
			dencount = dencount+1;
			name = tokens{5}(comma:end-12);
			
			if strcmp(tokens{2}, '0')
				zer(zercount) = j;
				zercount = zercount+1;
			end
		end
	end
	
	initial = newDendrite(1);
	newDendrite = newDendrite-initial+1;
	zer = zer - initial+1; % Where are the zeros?
	
	numberOfZeros = length(zer);
	newDendrite = [newDendrite (totalFiles+1)]; % Where are the new files?
		
	Numb = length(newDendrite)-1;
	state.imageProc.spineData.numberOfNeurons = Numb;
	updateGUIByGlobal('state.imageProc.spineData.numberOfNeurons');
	
	
	counter = 1;
	dendriteCounter=1;
	fileCounter=2;
	for i = 1:length(text)
		index = findstr('Dendrite	Spine Density(2D)' , text{i});
		if ~isempty(index)
			for j = 1:totalFiles
				tokens = tokenize(text{i+j});
				if str2num(tokens{1}) <= (newDendrite(fileCounter)-1)
					den2d{dendriteCounter}(counter) = str2num(tokens{2});
					counter = counter + 1;
				else
					counter = 1;
					dendriteCounter = dendriteCounter+1;
					fileCounter = fileCounter+1;
					den2d{dendriteCounter}(counter) = str2num(tokens{2});
					counter = counter+1;
				end
			end
		end	
	end
	
	for i = 1:length(den2d)
		values = den2d{i}(1:end);
		mean2d(i) = mean(values);
	end
	den2derr = std(mean2d);
	mean2d = mean(mean2d);
	
	counter = 1;
	dendriteCounter=1;
	fileCounter=2;
	for i = 1:length(text)
		index = findstr('Dendrite	Spine Density(3D)' , text{i});
		if ~isempty(index)
			for j = 1:totalFiles
				tokens = tokenize(text{i+j});
				if str2num(tokens{1}) <= (newDendrite(fileCounter)-1)
					den3d{dendriteCounter}(counter) = str2num(tokens{2});
					counter = counter + 1;
				else
					counter = 1;
					dendriteCounter = dendriteCounter+1;
					fileCounter = fileCounter+1;
					den3d{dendriteCounter}(counter) = str2num(tokens{2});
					counter = counter+1;
				end
			end
		end	
	end
	
	for i = 1:length(den3d)
		values = den3d{i}(1:end);
		mean3d(i) = mean(values);
	end
	den3derr = std(mean3d);
	mean3d = mean(mean3d);
	
	totalFiles = totalFiles - numberOfZeros;
	% Redefine the number of files subtracting the zero length contributiopns.
	% Recaclualte where the new dendritesa begin acconting for the zeros
	for i = 2:length(newDendrite)
		for j = 1:length(zer)
			if zer(j) <= newDendrite(i) & zer(j) >= newDendrite(i-1) 
				newDendrite(i) = newDendrite(i) - 1;
			end
		end
	end
	newDendrite(end) = 	totalFiles+1;
	
	counter = 1;
	dendriteCounter=1;
	fileCounter=2;
	
	for i = 1:length(text)
		index = findstr('Dendrite	Mean length' , text{i});
		if ~isempty(index)
			for j = 1:totalFiles
				tokens = tokenize(text{i+j});
				if str2num(tokens{1}) <= (newDendrite(fileCounter)-1)
					len2d{dendriteCounter}(counter) = str2num(tokens{2});
					counter = counter + 1;
				else
					counter = 1;
					dendriteCounter = dendriteCounter+1;
					fileCounter = fileCounter+1;
					len2d{dendriteCounter}(counter) = str2num(tokens{2});
					counter = counter+1;
				end
			end
		end	
	end

	for i = 1:length(len2d)
		values = len2d{i}(1:end);
		meanlen2d(i) = mean(values);
	end
	len2derr = std(meanlen2d);
	meanlen2d = mean(meanlen2d);
	
	counter = 1;
	dendriteCounter=1;
	fileCounter=2;
	for i = 1:length(text)
		index = findstr('Dendrite	Mean volume' , text{i});
		if ~isempty(index)
			for j = 1:totalFiles
				tokens = tokenize(text{i+j});
				if str2num(tokens{1}) <= (newDendrite(fileCounter)-1)
					vol2d{dendriteCounter}(counter) = str2num(tokens{2});
					counter = counter + 1;
				else
					counter = 1;
					dendriteCounter = dendriteCounter+1;
					fileCounter = fileCounter+1;
					vol2d{dendriteCounter}(counter) = str2num(tokens{2});
					counter = counter+1;
				end
			end
		end	
	end
	
	for i = 1:length(vol2d)
		values = vol2d{i}(1:end);
		meanvol2d(i) = mean(values);
	end
	vol2derr = std(meanvol2d);
	meanvol2d = mean(meanvol2d);
	
	meanSpineLength = meanlen2d;
	meanSpineLenError = len2derr;
	meanSpineVolume = meanvol2d;
	meanSpineVolError = vol2derr;
	densityNumber = [mean2d mean3d];
	densityError = [den2derr den3derr];
end


