function plotLines
	global FWHMdata offsetWave roi_1x roi_1y roi_2x roi_2y
	plot(FWHMdata);
	append(offsetWave);
	appendxy(roi_1x, roi_1y);
	appendxy(roi_2x, roi_2y);