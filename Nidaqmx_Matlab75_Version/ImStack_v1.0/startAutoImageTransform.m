function startAutoImageTransform
global state gh


if state.imageProc.autoSelect == 1
	currentfile = 1;
	lastfile = length(state.imageProc.names);
else
	currentfile = state.imageProc.auto.startLoadImage;
	lastfile = state.imageProc.auto.endLoadImage;
end

% Making and loading the input file
h = waitbar(0,'Starting Batch Image Processing....', 'Name', 'Batch Processing Status', ...
	'Position', [164.250000000000   693.750000000000   270.000000000000   056.250000000000], ...
	'Pointer', 'watch');

counter = 0;
avgImage=[];
waitTillAllImgAreDone=0;
for k = currentfile:lastfile
	counter = counter + 1;
	% determine the file to open
	if state.imageProc.autoSelect == 0
		state.imageProc.auto.startLoadImage = k;
		updateGUIByGlobal('state.imageProc.auto.startLoadImage');
		path = get(gh.autotransformGUI.pathLoadName, 'String');
		base = get(gh.autotransformGUI.baseLoadName, 'String');	
		
		filetype = get(gh.autotransformGUI.fileType, 'String');	
		value = get(gh.autotransformGUI.fileType, 'Value');
		filetype = filetype{value};
		
		if state.imageProc.auto.appendZeros == 1
			if k < 10
				base = [base '0' num2str(k)];
			else
				base = [base num2str(k)];
			end
		else
			base = [base num2str(k)];
		end
		
		fullname = [path base filetype];
	else
		fullname = [state.imageProc.pathName state.imageProc.names{k}];
		filetype = get(gh.autotransformGUI.fileType, 'String');	
		value = get(gh.autotransformGUI.fileType, 'Value');
		filetype = filetype{value};
	end
	% determine the file type to do
	if state.imageProc.autoSelect
		waitbar(counter/(lastfile-currentfile + 1),h, ['Doing file ' state.imageProc.names{k}]);
	else
		waitbar(counter/(lastfile-currentfile + 1),h, ['Doing file ' num2str(k) ' of ' num2str(lastfile)]);
	end
	
	switch filetype
	case '.tif'
		tempimage = opentif(fullname);
		header = readImageHeaderTif(fullname); % Read the image header
	case '.cfd'
		tempimage = openACFDSpine(fullname,state.imageProc.cfd.numberofChannels)
	case '.CFD'
		tempimage = openACFDSpine(fullname,state.imageProc.cfd.numberofChannels)
	otherwise
		disp('File Not Available');
	end
	
	
	% Doing the transform
	
	transform = get(gh.autotransformGUI.transform, 'String');
	transformvalue = get(gh.autotransformGUI.transform, 'Value');
	transform = transform{transformvalue};
	collected=0;
	averaged=0;
	projected=0;
	
	switch transform
	case 'Collect Across Images'
		try
			if counter > 1 
				avgImage = cat(3,avgImage,double(tempimage));
			else
				avgImage = double(tempimage);
			end
		catch
			disp('Images not the same size. Cant collect');
            disp(lasterr);
		end
		cl=class(tempimage);
		extracted = 1;
		waitTillAllImgAreDone=1;	% Need to load when done all images being loaded
		collected=1;
	case 'Project Across Images'
		try
			if counter > 1 
				avgImage = collapse(cat(3,avgImage,double(tempimage)),'XY');
			else
				avgImage = double(tempimage);
			end
		catch
			disp('Images not the same size. Cant project.');
            disp(lasterr);
		end
		cl=class(tempimage);
		extracted = 1;
		waitTillAllImgAreDone=1;	% Need to load when done all images being loaded	
		projected=1;

	case 'Average Across Images'
		try
			if counter > 1 
				avgImage = avgImage + double(tempimage);
			else
				avgImage = double(tempimage);
			end
		catch
			disp('Images not the same size. Cant average');
            disp(lasterr);
		end
		cl=class(tempimage);
		extracted = 1;
		waitTillAllImgAreDone=1;	% Need to load when done all images being loaded
		averaged=1;
	case 'Collect and Filter Across Images'
		try
			if counter > 1 
				avgImage = cat(3,avgImage,medfilt2(double(tempimage),[3 3]));
			else
				avgImage =  medfilt2(double(tempimage),[3 3]);
			end
		catch
			disp('Images not the same size. Cant collect');
            disp(lasterr);
		end
		cl=class(tempimage);
		extracted = 1;
		waitTillAllImgAreDone=1;	% Need to load when done all images being loaded
		collected=1;
	case 'Max Projection (XY)'
		tempimage = collapse(tempimage, 'XY');
		extracted = 0;
	case 'Max Projection (XZ)'
		tempimage = collapse(tempimage, 'XZ');
		extracted = 0;
	case 'Max Projection (YZ)'
		tempimage = collapse(tempimage, 'YZ');
		extracted = 0;
	case 'Extract and Project XY'
		tempimage = extractData(tempimage, state.imageProc.cfd.numberofChannels);
		for i = 1:state.imageProc.cfd.numberofChannels
			tempimage{i} = collapse(tempimage{i}, 'XY');
			
			if state.imageProc.autoSelect == 0
				savepath = get(gh.autotransformGUI.pathSaveName, 'String');
				savebase = get(gh.autotransformGUI.baseSaveName, 'String');
				currentsavefile = get(gh.autotransformGUI.startSaveImage, 'String');
				
				if state.imageProc.auto.appendZerosSave == 1
					if state.imageProc.auto.startSaveImage < 10
						fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) 'max00' currentsavefile ];
					else
						fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) 'max0' currentsavefile ];
					end
				else
					fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) 'max' currentsavefile ];
				end
			else
				savepath = get(gh.autotransformGUI.pathSaveName, 'String');
				[path,name,ext] = fileparts(state.imageProc.names{k});
				fullname = [savepath num2str(i) name];
			end
			
			arrayToTiff(tempimage{i}, fullname, header);
			
			% Should we load the file?
			
			if state.imageProc.auto.loadAsTransform == 1
				[path, name, ext] = fileparts(fullname);
				loadImageFromName([path '\' name '.tif'], path);
			end
		end
		if state.imageProc.autoSelect == 0
			state.imageProc.auto.startSaveImage = str2num(currentsavefile) + 1;
			updateGUIByGlobal('state.imageProc.auto.startSaveImage');
		end
		
	case 'Extract'
		tempimage = extractData(tempimage, state.imageProc.cfd.numberofChannels);
		
		for i = 1:state.imageProc.cfd.numberofChannels
			tempimage{i} = collapse(tempimage{i}, 'XY');
			
			if state.imageProc.autoSelect == 0
				savepath = get(gh.autotransformGUI.pathSaveName, 'String');
				savebase = get(gh.autotransformGUI.baseSaveName, 'String');
				currentsavefile = get(gh.autotransformGUI.startSaveImage, 'String');
				
				if state.imageProc.auto.appendZerosSave == 1
					if state.imageProc.auto.startSaveImage < 10
						fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) '00' currentsavefile ];
					else
						fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) '0' currentsavefile ];
					end
				else
					fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) currentsavefile ];
				end
			else
				savepath = get(gh.autotransformGUI.pathSaveName, 'String');
				[path,name,ext] = fileparts(state.imageProc.names{k});
				fullname = [savepath num2str(i) name];
			end
			
			arrayToTiff(tempimage{i}, fullname, header);
			
			% Should we load the file?
			
			if state.imageProc.auto.loadAsTransform == 1
				[path, name, ext] = fileparts(fullname);
				loadImageFromName([path '\' name '.tif'], path);
			end
		end
		if state.imageProc.autoSelect == 0
			state.imageProc.auto.startSaveImage = str2num(currentsavefile) + 1;
			updateGUIByGlobal('state.imageProc.auto.startSaveImage');
		end
		extracted = 1;
	case 'Filter and Project XYZ'
		for j = 1:size(tempimage,3)
			tempimage(:,:,j) = medfilt2(tempimage(:,:,j),[3 3]);
		end
		tempimage1{1} = collapse(tempimage, 'XY');
		tempimage1{2} = collapse(tempimage, 'XZ');
		tempimage1{3} = collapse(tempimage, 'YZ');
		tempimage = tempimage1;
		
		for i = 1:3
			switch i 
			case 1
				dir = 'maxXY';
			case 2
				dir = 'maxXZ';
			case 3
				dir = 'maxYZ';
			end
			
			if state.imageProc.autoSelect == 0
				savepath = get(gh.autotransformGUI.pathSaveName, 'String');
				savebase = get(gh.autotransformGUI.baseSaveName, 'String');
				currentsavefile = get(gh.autotransformGUI.startSaveImage, 'String');
				
				if state.imageProc.auto.appendZerosSave == 1
					if state.imageProc.auto.startSaveImage < 10
						fullname = [savepath savebase(1:(end-1)) dir '00' currentsavefile ];
					else
						fullname = [savepath savebase(1:(end-1)) dir '0' currentsavefile ];
					end
				else
					fullname = [savepath savebase(1:(end-1))  dir currentsavefile ];
				end
			else
				savepath = get(gh.autotransformGUI.pathSaveName, 'String');
				[path,name,ext] = fileparts(state.imageProc.names{k});
				fullname = [savepath dir name];
			end
			
			arrayToTiff(tempimage{i}, fullname, header);
			
			% Should we load the file?
			
			if state.imageProc.auto.loadAsTransform == 1
				[path, name, ext] = fileparts(fullname);
				loadImageFromName([path '\' name '.tif'], path);
			end
		end
		if state.imageProc.autoSelect == 0
			state.imageProc.auto.startSaveImage = str2num(currentsavefile) + 1;
			updateGUIByGlobal('state.imageProc.auto.startSaveImage');
		end
		extracted = 1;
	case 'Extract and Filter'
		tempimage = extractData(tempimage, state.imageProc.cfd.numberofChannels);
		for i = 1:state.imageProc.cfd.numberofChannels
			for j = 1:size(tempimage{i},3)
				tempimage{i}(:,:,j) = medfilt2(tempimage{i}(:,:,j),[3 3]);
			end
			if state.imageProc.autoSelect == 0
				savepath = get(gh.autotransformGUI.pathSaveName, 'String');
				savebase = get(gh.autotransformGUI.baseSaveName, 'String');
				currentsavefile = get(gh.autotransformGUI.startSaveImage, 'String');
				
				if state.imageProc.auto.appendZerosSave == 1
					if state.imageProc.auto.startSaveImage < 10
						fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) 'f00' currentsavefile ];
					else
						fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) 'f0' currentsavefile ];
					end
				else
					fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) 'f' currentsavefile ];
				end
			else
				savepath = get(gh.autotransformGUI.pathSaveName, 'String');
				[path,name,ext] = fileparts(state.imageProc.names{k});
				fullname = [savepath num2str(i) name];
			end
			
			arrayToTiff(tempimage{i}, fullname, header);
			
			% Should we load the file?
			
			if state.imageProc.auto.loadAsTransform == 1
				[path, name, ext] = fileparts(fullname);
				loadImageFromName([path '\' name '.tif'], path);
			end
		end
		if state.imageProc.autoSelect == 0
			state.imageProc.auto.startSaveImage = str2num(currentsavefile) + 1;
			updateGUIByGlobal('state.imageProc.auto.startSaveImage');
		end
		
		extracted = 1;
	case 'Extract, Filter, & Project'
		tempimage = extractData(tempimage, state.imageProc.cfd.numberofChannels);
		for i = 1:state.imageProc.cfd.numberofChannels
			for j = 1:size(tempimage{i},3)
				tempimage{i}(:,:,j) = medfilt2(tempimage{i}(:,:,j),[3 3]);
			end
			tempimage{i} = collapse(tempimage{i}, 'XY');
			if state.imageProc.autoSelect == 0
				savepath = get(gh.autotransformGUI.pathSaveName, 'String');
				savebase = get(gh.autotransformGUI.baseSaveName, 'String');
				currentsavefile = get(gh.autotransformGUI.startSaveImage, 'String');
				
				if state.imageProc.auto.appendZerosSave == 1
					if state.imageProc.auto.startSaveImage < 10
						fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) 'fmax00' currentsavefile ];
					else
						fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) 'fmax0' currentsavefile ];
					end
				else
					fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) 'fmax' currentsavefile ];
				end
			else
				savepath = get(gh.autotransformGUI.pathSaveName, 'String');
				[path,name,ext] = fileparts(state.imageProc.names{k});
				fullname = [savepath num2str(i) name];
			end
			
			arrayToTiff(tempimage{i}, fullname, header);
			
			% Should we load the file?
			
			if state.imageProc.auto.loadAsTransform == 1
				[path, name, ext] = fileparts(fullname);
				loadImageFromName([path '\' name '.tif'], path);
			end
		end
		if state.imageProc.autoSelect == 0
			state.imageProc.auto.startSaveImage = str2num(currentsavefile) + 1;
			updateGUIByGlobal('state.imageProc.auto.startSaveImage');
		end
		
		extracted = 1;
	case 'Re-Write Header'
		%header = changeHeaderFields(header);
		extracted = 0;
	case 'Filter'
		for i = 1:size(tempimage,3)
			tempimage(:,:,i) = medfilt2(tempimage(:,:,i),[3 3]);
		end
		extracted = 0;
	case 'Filter and Project XY'
		for i = 1:size(tempimage,3)
			tempimage(:,:,i) = medfilt2(tempimage(:,:,i),[3 3]);
		end
		tempimage = collapse(tempimage, 'XY');
		extracted = 0;
	case 'Average'
		tempimage = genericBin(tempimage, 1,1,size(tempimage,3));
		extracted = 0;
	case 'Extract and Average'
		tempimage = extractData(tempimage, state.imageProc.cfd.numberofChannels);
		for i = 1:state.imageProc.cfd.numberofChannels
			tempimage{i} = genericBin(tempimage{i}, 1,1,size(tempimage{i},3));
			if state.imageProc.autoSelect == 0
				savepath = get(gh.autotransformGUI.pathSaveName, 'String');
				savebase = get(gh.autotransformGUI.baseSaveName, 'String');
				currentsavefile = get(gh.autotransformGUI.startSaveImage, 'String');
				
				if state.imageProc.auto.appendZerosSave == 1
					if state.imageProc.auto.startSaveImage < 10
						fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) 'avg00' currentsavefile ];
					else
						fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) 'avg0' currentsavefile ];
					end
				else
					fullname = [savepath savebase(1:(end-1)) 'c' num2str(i) 'avg' currentsavefile ];
				end
			else
				savepath = get(gh.autotransformGUI.pathSaveName, 'String');
				[path,name,ext] = fileparts(state.imageProc.names{k});
				fullname = [savepath num2str(i) name];
			end
			
			arrayToTiff(tempimage{i}, fullname, header);
			
			% Should we load the file?
			
			if state.imageProc.auto.loadAsTransform == 1
				[path, name, ext] = fileparts(fullname);
				loadImageFromName([path '\' name '.tif'], path);
			end
		end
		if state.imageProc.autoSelect == 0
			state.imageProc.auto.startSaveImage = str2num(currentsavefile) + 1;
			updateGUIByGlobal('state.imageProc.auto.startSaveImage');
		end
		extracted = 1;
		
	case 'Extract, Average, Ratio, & Threshold'
		tempimage = extractData(tempimage, state.imageProc.cfd.numberofChannels);
		for i = 1:state.imageProc.cfd.numberofChannels
			for j = 1:size(tempimage{i},3)
				tempimage{i}(:,:,j) = medfilt2(tempimage{i}(:,:,j),[3 3]);
			end
		end
		%for i = 1:state.imageProc.cfd.numberofChannels
		%	tempimage{i} = genericBin(tempimage{i}, 1,1,size(tempimage{i},3));
		%end
		
		%binary = tempimage{1} > state.imageProc.parsing.threshold;
		%tempimage{1} = double(tempimage{1}).*double(binary));
		%tempimage{2} =(double(tempimage{2}).*double(binary));
		tempimage =double(tempimage{1})./double(tempimage{2});
		[r,c] = find(tempimage == inf);
		for j = 1:length(r)
			tempimage(r(j),c(j)) = state.imageProc.parsing.threshold+1;
		end
		[r,c] = find(tempimage == NaN);
		for j = 1:length(r)
			tempimage(r(j),c(j)) = 0;
		end
		
		extracted = 0;
	otherwise
		disp('No transform selected');
	end
	
	% Saving the new array as a tif
	
	
	if ~extracted 
		if state.imageProc.autoSelect == 0
			savepath = get(gh.autotransformGUI.pathSaveName, 'String');
			savebase = get(gh.autotransformGUI.baseSaveName, 'String');
			currentsavefile = get(gh.autotransformGUI.startSaveImage, 'String');
			
			if state.imageProc.auto.appendZerosSave == 1
				if state.imageProc.auto.startSaveImage < 10
					fullname = [savepath savebase '0' currentsavefile];
				else
					fullname = [savepath savebase currentsavefile];
				end
			else
				fullname = [savepath savebase currentsavefile];
			end
		else
			savepath = get(gh.autotransformGUI.pathSaveName, 'String');
			[path,name,ext] = fileparts(state.imageProc.names{k});
			fullname = [savepath name];
		end
		
		arrayToTiff(tempimage, fullname, header);
		% Should we load the file?
		
		if state.imageProc.auto.loadAsTransform == 1
			[path, name, ext] = fileparts(fullname);
			loadImageFromName([path '\' name '.tif'], path);
		else
		end
		if state.imageProc.autoSelect == 0
			state.imageProc.auto.startSaveImage = str2num(currentsavefile) + 1;
			updateGUIByGlobal('state.imageProc.auto.startSaveImage');
		end
	end
end

%now doen loading all the files so can do some global stuff...

if waitTillAllImgAreDone
	if state.imageProc.autoSelect == 0
		savepath = get(gh.autotransformGUI.pathSaveName, 'String');
		savebase = get(gh.autotransformGUI.baseSaveName, 'String');
		currentsavefile = get(gh.autotransformGUI.startSaveImage, 'String');
		
		if state.imageProc.auto.appendZerosSave == 1
			if state.imageProc.auto.startSaveImage < 10
				fullname = [savepath savebase '0' currentsavefile];
			else
				fullname = [savepath savebase currentsavefile];
			end
		else
			fullname = [savepath savebase currentsavefile];
		end
	else
		savepath = get(gh.autotransformGUI.pathSaveName, 'String');
		[path,name,ext] = fileparts(state.imageProc.names{k});
		fullname = [savepath name];
	end
	if collected | projected
		%bdont do anything.
	end
	if 	averaged
		avgImage = 1/(counter-1)*avgImage;	%do avg.
	end
	arrayToTiff(eval([cl '(avgImage);']), fullname, header);
	% Should we load the file?
	if state.imageProc.auto.loadAsTransform == 1
		[path, name, ext] = fileparts(fullname);
		loadImageFromName([path '\' name '.tif'], path);
	end
end
close(h);
	
	
	
