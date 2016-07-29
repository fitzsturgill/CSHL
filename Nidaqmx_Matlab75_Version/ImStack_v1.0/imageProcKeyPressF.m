function imageProcKeyPressF

val = double(get(gcbo,'CurrentCharacter'));

try
	switch val
		case 6 % ctrl + f
			applyMedianFilter;
		case 14 % ctrl + n
			loadImage;
		case 13 %ctrl + m
			loadProjection;
		case 4 %ctrl + d
			closeAllImages;
		case 1 % ctrl + a
			computeTimeSeriesMean;
		case 24 % ctrl + x
			computeTimeSeriesMeanV;
		case 26 % ctrl + z
			loadROI;
		case 29
			global state gh

			pathName = get(gh.fileCounterGUI.pathName, 'String');
			baseName = get(gh.fileCounterGUI.baseName, 'String');

	
			value = get(gh.fileCounterGUI.filetype, 'Value');
			ext = get(gh.fileCounterGUI.filetype, 'String');
			ext = ext{value};
			acquisitionNumber = get(gh.fileCounterGUI.acquisitionCounter, 'String');

			if state.imageProc.appendZeros == 1
				if str2num(acquisitionNumber) < 10
					baseName = [baseName '0'];
				else
				end
			end

			filename = [pathName baseName num2str(acquisitionNumber) ext];
			if strcmp(ext, '.tif')
				loadImageFromName(filename,pathName);
			elseif strcmp(ext, '.jpg')
				state.imageProc.currentjpeg = imread(filename);
				state.imageProc.currentjpeg = rgb2gray(state.imageProc.currentjpeg);
				loadImageFromArray('state.imageProc.currentjpeg');
    		elseif strcmp(ext, '.CFD')
				seeGUI('gh.cfdGUI.figure1');
				set(gh.fileCounterGUI.filetype, 'Value', 3);
    		    openACFD(filename);
				hideGUI('gh.cfdGUI.figure1');
    		elseif strcmp(ext, '.cfd')
				set(gh.fileCounterGUI.filetype, 'Value', 3);
  		      seeGUI('gh.cfdGUI.figure1');
   		     openACFD(filename);
			hideGUI('gh.cfdGUI.figure1');
			end
			state.imageProc.acquisitionNumber = state.imageProc.acquisitionNumber+1;
			updateGUIByGlobal('state.imageProc.acquisitionNumber');


	otherwise
	end
end
