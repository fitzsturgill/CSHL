function recordMovie(type)
global gh state

if nargin < 1
	
	[fname pname] = uiputfile('*.avi', 'Save movie as...');
	
	[path,name,ext] = fileparts([pname fname]);
	
	if pname < 1
		return
	end
	
	filename = [pname fname];
	state.imageProc.currentMovie = [];
	
	value = get(gh.movieGUI.fileName, 'Value');
	name = get(gh.movieGUI.fileName, 'String');
	if iscell(name)
		name = name{value};
	else
	end
	
	array = state.imageProc.cell.currentImage{value}(:,:,state.imageProc.movieStart:state.imageProc.movieEnd);
	state.imageProc.movieFigure = figure('NumberTitle', 'off', 'Name', ['Movie of ' name ' from ' num2str(state.imageProc.movieStart) ...
			' to ' num2str(state.imageProc.movieEnd) '.'], 'Position', [360 360 size(array,1) size(array,2)],'DoubleBuffer', 'on');
	set(gca,'nextplot','replacechildren', 'YDir', 'Reverse', 'Position', [0 0 1 1] , 'XTickLabel', [], 'YTickLabel', [], ...
		'XTickLabelMode', 'manual', 'YTickLabelMode', 'manual', 'visible', 'off');
	
	frames = size(array,3);
	array = double(array);
	
	mov = avifile(filename, 'fps', state.imageProc.fps, 'compression', 'indeo5')
	
	% Record the movie
	for j = 1:frames 
		h = image('CData', array(:,:,j), 'CDataMapping', 'Scaled');
		set(gca, 'CLim', [state.imageProc.cell.lowPixelValue{value} state.imageProc.cell.highPixelValue{value}]);
		axis image;
		colormap(gray);
		F(j) = getframe(gca);
		mov = addframe(mov,F(j));
	end
	state.imageProc.currentMovie = F;
	mov = close(mov);
	
else
	if strcmp(type, 'matlab')
		
		value = get(gh.movieGUI.fileName, 'Value');
		name = get(gh.movieGUI.fileName, 'String');
		if iscell(name)
			name = name{value};
		else
		end
	
		array = state.imageProc.cell.currentImage{value}(:,:,state.imageProc.movieStart:state.imageProc.movieEnd);
		state.imageProc.movieFigure = figure('NumberTitle', 'off', 'Name', ['Movie of ' name ' from ' num2str(state.imageProc.movieStart) ...
				' to ' num2str(state.imageProc.movieEnd) '.'], 'Position', [360 360 size(array,1) size(array,2)],'DoubleBuffer', 'on');
		set(gca,'nextplot','replacechildren', 'YDir', 'Reverse', 'Position', [0 0 1 1] , 'XTickLabel', [], 'YTickLabel', [], ...
			'XTickLabelMode', 'manual', 'YTickLabelMode', 'manual', 'visible', 'off');
		
		frames = size(array,3);
		array = double(array);
				
		% Record the movie
		for j = 1:frames 
			h = image('CData', array(:,:,j), 'CDataMapping', 'Scaled');
			set(gca, 'CLim', [state.imageProc.cell.lowPixelValue{value} state.imageProc.cell.highPixelValue{value}]);
			axis image;
			colormap(gray);
			F(j) = getframe(gca);
		end
		state.imageProc.currentMovie = F;
	else 
		disp('Error: only valid paraemter is matlab.');
	end
end

