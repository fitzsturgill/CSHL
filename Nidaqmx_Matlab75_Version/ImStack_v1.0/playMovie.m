function playMovie
global state 


% Play current Movie
try
	if state.imageProc.movieFigure < 1
		state.imageProc.movieFigure = figure('NumberTitle', 'off', 'Name', ['Movie of ' name ' from ' num2str(state.imageProc.movieStart) ...
		' to ' num2str(state.imageProc.movieEnd) '.'],'DoubleBuffer', 'on');
		movie(gca, state.imageProc.currentMovie,state.imageProc.repeat,state.imageProc.fps);
		axis image;
		colormap(gray);
	else
		set(state.imageProc.movieFigure, 'Visible', 'On');
		movie(gca,state.imageProc.currentMovie,state.imageProc.repeat,state.imageProc.fps);
		axis image;
		colormap(gray);
	end
catch 
end
