function currentIntensity = getCurrentIntensity(axis, image)
	
	% this function returns the value or inensisty at a point last cliked in the image.
	
	currentPoint = recordCurrentPoint(axis);
	CData = get(image, 'CData');
	currentIntensity = CData(currentPoint(2), currentPoint(1));
	