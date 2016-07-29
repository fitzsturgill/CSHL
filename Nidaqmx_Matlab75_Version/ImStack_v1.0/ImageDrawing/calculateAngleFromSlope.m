function angle = calculateAngleFromSlope(object)

type = get(object, 'Type');

try
	if strcmp(type, 'line')
		
		Ydata = get(object, 'YData');
		Xdata = get(object, 'XData');

		slope = (Ydata(2)-Ydata(1))/(Xdata(2)-Xdata(1));
		angleDegrees = 180*atan(slope)/pi;
		angleRads = atan(slope);
		angle = [angleDegrees angleRads];
	else
		angle = [];
		display('Object input must be a line object.');
	end
	
catch
	angle = [];
	display('Object input must be a line object.');
end
