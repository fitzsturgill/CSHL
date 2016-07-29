function loadImage
global state gh

% This function will call a ui to open a tif file and set it up for reading.
% It will aslo update the glboals associated with that particular image.

% Call load image ui
[fname, pname] = uigetfile({'*.tif;*.jpg;*.cfd;*.CFD'}, 'Choose image to load');
if fname < 1
	return
else

	set(gh.fileCounterGUI.baseName, 'String', fname);
	set(gh.fileCounterGUI.pathName, 'String', pname);
	[path,name,ext] = fileparts([pname fname]);
	
	switch ext 
		case '.jpg'
			fullname=[path '\' name ext];
			state.imageProc.currentjpeg = [];
			state.imageProc.currentjpeg = imread(fullname);
			try
				state.imageProc.currentjpeg = rgb2gray(state.imageProc.currentjpeg);
			catch
			end	
			loadImageFromArray('state.imageProc.currentjpeg');
			set(gh.fileCounterGUI.filetype, 'Value', 2);
            try
			    cd(pname);
            end
		case '.CFD'
			set(gh.fileCounterGUI.filetype, 'Value', 3);
			state.imageProc.colorMap=0;
  			loadImageFromName([pname fname], pname);
            try
			    cd(pname);
            end
			
		case '.cfd'
			set(gh.fileCounterGUI.filetype, 'Value', 3);
			state.imageProc.colorMap=0;
      	  	openACFD([pname fname]);
		
		case '.tif'
			set(gh.fileCounterGUI.filetype, 'Value', 1);
  			loadImageFromName([pname fname], pname);
            try
			    cd(pname);
            end
            
			
	end
	
end

	
