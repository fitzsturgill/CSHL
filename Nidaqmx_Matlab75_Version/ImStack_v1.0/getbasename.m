function basename = getbasename(filename, pathname)
global gh state

lengthfile = length(filename);
lengthpath = length(pathname);
basename = filename((lengthpath+1):lengthfile);
filesize = length(basename);
periods=0;

if state.imageProc.appendZeros == 1 | state.imageProc.auto.appendZeros == 1
	
	for i=1:filesize
		if strcmp(basename(i), '.')
			basename = basename(1:(i-3));
			return
		end
	end

else
	
	for i=1:filesize
		if strcmp(basename(i), '.')
			basename = basename(1:(i-1));
			return
		end
	end
end


