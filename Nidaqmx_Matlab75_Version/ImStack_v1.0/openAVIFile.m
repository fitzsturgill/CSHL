function openAVIFile
global gh state

[fname pname] = uigetfile('*.avi', 'Select movie to open...');

if pname < 1
	return
else
		state.imageProc.movieAviFigure = figure('NumberTitle', 'off', 'Name', ['AVI of ' fname]);
		state.imageproc.currentMovie = aviread([pname fname]);
		axis image;
		set(gca,'nextplot','replacechildren', 'YDir', 'Reverse', 'Position', [0 0 1 1] , 'XTickLabel', [], 'YTickLabel', [], ...
			'XTickLabelMode', 'manual', 'YTickLabelMode', 'manual', 'visible', 'off');
		colormap(gray);
		movie(gca, state.imageProc.currentMovie);
end

	
	