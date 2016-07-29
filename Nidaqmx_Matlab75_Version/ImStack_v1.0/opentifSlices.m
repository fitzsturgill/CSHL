function Aout = opentif(filename)
global state

% Open a tif file and store its contents as array Aout.
% filename is the file name with extension.  
% 
% Written By: Thomas Pologruto
% Cold Spring Harbor Labs
% February 1, 2001

h = waitbar(0,'Opening Tif image...', 'Name', 'Open TIF Image', 'Pointer', 'watch');
[path,name,ext] = fileparts(filename);
cd(path);
sliceNames=dir([name(1:end-6) '*.tif']);
frames=size(sliceNames,1);
state.imageProc.colorMap = 0;
disp('Loading by slices');

	for i = 1:frames
		waitbar(i/frames,h, ['Loading Frame Number ' num2str(i)]);
		filename=sliceNames(i).name;
		disp(filename)
		try
			Aout(:,:,i) = imread(filename, 1);
		catch
			if i == 1  %If it is mismatch error on first frame, it is a rgb tif file.
				try
					Aout = imread(filename);
					state.imageProc.colorMap = 1;
					waitbar(1,h, 'Done');
					close(h);
					return
				catch
					beep;
					waitbar(1,h, 'Done');
					close(h);
					Aout = 'File Not Found';
					disp(['File ' filename ' Not Found; Check Name or Path.']);
				end
			end
		end	
	end

	waitbar(1,h, 'Done');
	close(h);



