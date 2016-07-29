function array = openACFDSpine(filename,channels)
global state gh head_info

head_info = CF4header(filename);
totalimages=head_info.n_images;

h = waitbar(0,'Opening CFD image...', 'Name', 'Open CFD Image', 'Pointer', 'watch');

try
	for i=1:totalimages	
	  waitbar(i/totalimages,h, ['Loading Frame Number ' num2str(i)]);
 	  array(:,:,i) = CFDread2(filename,i, channels);
	end
catch
	  waitbar(1,h, 'Done');
end   

close(h);