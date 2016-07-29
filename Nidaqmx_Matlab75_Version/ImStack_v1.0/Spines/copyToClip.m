function copyToClip(handle)
global gh state

	type = get(handle, 'Type');

	switch type
	case 'figure'
		set(handle, 'visible', 'on','InvertHardCopy','off');
		print(handle, '-dbitmap');
	case 'axes'
		XLim = get(handle, 'XLim');
		YLim = get(handle, 'YLim');
		x = XLim(2) - XLim(1);
		y = YLim(2) - YLim(1);
		f = figure('Position', [360 360 x y], 'menubar', 'none', 'numberTitle', 'off');
		colormap(gray);
		axis = copyobj(handle, f);
		set(axis, 'units', 'normalized', 'position', [0 0 1 1]);
		print(f, '-dbitmap');
		close(f);
	case 'image'
		parent=get(handle,'Parent');
		Ydir=get(parent,'YDir');
		Xdir=get(parent,'XDir');
		
		XLim = get(handle, 'XData');
		YLim = get(handle, 'YData');
		x = XLim(2) - XLim(1);
		y = YLim(2) - YLim(1);
		f = figure('Position', [360 360 x y], 'menubar', 'none', 'numberTitle', 'off');
		colormap(gray);
		image = copyobj(handle, gca);
		set(gca, 'units', 'normalized', 'position', [0 0 1 1], 'XLim', XLim', 'YLim', YLim,'Ydir',Ydir,...
			'Xdir',Xdir);
		print(f, '-dbitmap');
		close(f);
	otherwise
		disp('Handle invalid');
	end
