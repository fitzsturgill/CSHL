function saveMovieAsAVI(movie, fps)
global gh state

[fname pname] = uiputfile('*.avi', 'Save movie as...');

[path,name,ext] = fileparts([pname fname]);

if pname < 1
	return
else
	movie2avi(movie, [pname fname '.avi'], 'FPS', fps);
end
