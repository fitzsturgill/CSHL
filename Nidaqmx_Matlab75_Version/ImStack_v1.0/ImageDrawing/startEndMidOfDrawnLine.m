function [startEndMid] = startEndMidOfDrawnLine(axis, line)

% parameters pertaining to axis

position = get(axis, 'Position');
XLimits = get(axis, 'XLim');
YLimits = get(axis, 'YLim');
Ydirection = get(axis, 'YDir');

minValueX = position(1,1);
maxValueX = position(1,1) + position(1,3);
minValueY = position(1,2);
maxValueY = position(1,2) + position(1,4);
xLength = XLimits(1,2) - XLimits(1,1);

if strcmp(Ydirection, 'reverse')	
	yLength = YLimits(1,1) - YLimits(1,2); % Y axis is reverse
else
	yLength = YLimits(1,2) - YLimits(1,1); % Y axis is normal
	
end

% Parameters pertaining to the drawn line
YdataLine = get(line, 'YData');
XdataLine = get(line, 'XData');

% Conversion from axis to pixel coordinates is linear and follows the following formulae:

slopeX = xLength/(maxValueX-minValueX);
yintX = XLimits(1,1) - slopeX*minValueX;

slopeY = yLength/(maxValueY-minValueY);
yintY = YLimits(1,2) - slopeY*minValueY;


if XdataLine(1) <= maxValueX & XdataLine(2) <= maxValueX & YdataLine(1) <= maxValueY & ...
		YdataLine(2) <= maxValueY & XdataLine(1) >= minValueX & XdataLine(2) >= minValueX ...
		& YdataLine(1) >= minValueY & YdataLine(2) >= minValueY
	
	midpointX = .5*(XdataLine(1) + XdataLine(2));
	midpointX = round(slopeX*midpointX + yintX);
	
	midpointY = .5*(YdataLine(1) + YdataLine(2));
	midpointY = round(slopeY*midpointY + yintY);
	
	beginningPoint = [round(slopeX*XdataLine(1) + yintX) round(slopeY*YdataLine(1) + yintY)];
	endingPoint = [round(slopeX*XdataLine(2) + yintX) round(slopeY*YdataLine(2) + yintY)];
	midPoint = [midpointX, midpointY];
	
	startEndMid = [beginningPoint; endingPoint; midPoint];
else
	display('Line off image. Please Redraw.');
end
