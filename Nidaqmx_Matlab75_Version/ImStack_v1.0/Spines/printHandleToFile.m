function printHandleToFile(handle, filename)
global gh state

	type = get(handle, 'Type');

	switch type
	case 'figure'
		set(handle, 'Visible', 'on');
		set(handle, 'menubar', 'none', 'numberTitle', 'off', 'PaperPositionMode', 'auto', ...
			'InvertHardcopy','off');
		print(handle, '-dtiff', filename);
		set(handle, 'menubar', 'figure');
	case 'axes'
		XLim = get(handle, 'XLim');
		YLim = get(handle, 'YLim');
		x = XLim(2) - XLim(1);
		y = YLim(2) - YLim(1);
		f = figure('Position', [360 360 x y], 'menubar', 'none', 'numberTitle', 'off', 'PaperPositionMode', 'auto', ...
			'InvertHardcopy','off');
		colormap(gray);
		axis = copyobj(handle, f);
		set(axis, 'units', 'normalized', 'position', [0 0 1 1]);
		print(f, '-dtiff', filename);
		close(f);
	case 'image'
		XLim = get(handle, 'XData');
		YLim = get(handle, 'YData');
		x = XLim(2) - XLim(1);
		y = YLim(2) - YLim(1);
		f = figure('Position', [360 360 x y], 'menubar', 'none', 'numberTitle', 'off', 'PaperPositionMode', 'auto', ...
			'InvertHardcopy','off');
		colormap(gray);
		image = copyobj(handle, gca);
		set(gca, 'units', 'normalized', 'position', [0 0 1 1], 'XLim', XLim, 'YLim', YLim);
		print(f, '-dtiff', filename);
		close(f);
	otherwise
		disp('Handle invalid');
	end