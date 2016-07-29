function [annovaName, zerosRemovedAnova,  anDensity2d, anDensity3d, anLen, anVol]  = anSpineDensityVar(filename)
global gh state

text = textread(filename, '%s', 'headerlines', 1, 'delimiter', '\n');

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
anCount = 1;
name = 'test';
for j = start:endO
	tokens = tokenize(text{j});
	comma = findstr(',', tokens{5})+1;
	if strcmp(name, tokens{5}(comma:end-12)) 
		annovaName{anCount} = name;
		anCount = anCount+1;
		if strcmp(tokens{2}, '0')
			zer(zercount) = j;
			zercount = zercount+1;	
		end
	else
		newDendrite(dencount) = j;
		dencount = dencount+1;
		name = tokens{5}(comma:end-12);
		annovaName{anCount} = name;
		anCount = anCount+1;
		if strcmp(tokens{2}, '0')
			zer(zercount) = j;
			zercount = zercount+1;
		end
	end
end

initial = newDendrite(1);
newDendrite = newDendrite-initial+1;
zer = zer - initial+1; % Where are the zeros?

annovaName; % Theis is a cell array of the neuron names
zerosRemovedAnova = deleteCellContents(annovaName, zer); % cell array of neuron names with zero spine dendrites removed/

numberOfZeros = length(zer);
newDendrite = [newDendrite (totalFiles+1)]; % Where are the new files?

Numb = length(newDendrite)-1;
state.imageProc.spineData.numberOfNeurons = Numb;
updateGUIByGlobal('state.imageProc.spineData.numberOfNeurons');

for i = 1:length(text)
	index = findstr('Total no. of spines:'  , text{i});
	if ~isempty(index)
		tokens = tokenize(text{i});
		state.imageProc.spineData.numberOfSpines = str2num(tokens{5});
		updateGUIByGlobal('state.imageProc.spineData.numberOfSpines');
		break
		
	end	
end


counter = 1;
dendriteCounter=1;
anDens = 1;
fileCounter=2;
for i = 1:length(text)
	index = findstr('Dendrite	Spine Density(2D)' , text{i});
	if ~isempty(index)
		for j = 1:totalFiles
			tokens = tokenize(text{i+j});
			anDensity2d(anDens) = str2num(tokens{2});
			anDens = anDens+1;
		end
	end
end	


counter = 1;
dendriteCounter=1;
anDens = 1;
fileCounter=2;
for i = 1:length(text)
	index = findstr('Dendrite	Spine Density(3D)' , text{i});
	if ~isempty(index)
		for j = 1:totalFiles
			tokens = tokenize(text{i+j});
			anDensity3d(anDens) = str2num(tokens{2});
			anDens = anDens+1;
		end
	end
end	



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
anDens = 1;
fileCounter=2;
for i = 1:length(text)
	index = findstr('Dendrite	Mean length' , text{i});
	if ~isempty(index)
		for j = 1:totalFiles
			tokens = tokenize(text{i+j});
			anLen(anDens) = str2num(tokens{2});
			anDens = anDens+1;
		end
	end
end	



counter = 1;
dendriteCounter=1;
anDens = 1;
fileCounter=2;
for i = 1:length(text)
	index = findstr('Dendrite	Mean volume' , text{i});
	if ~isempty(index)
		for j = 1:totalFiles
			tokens = tokenize(text{i+j});
			anVol(anDens) = str2num(tokens{2});
			anDens = anDens+1;
		end
	end
end	



