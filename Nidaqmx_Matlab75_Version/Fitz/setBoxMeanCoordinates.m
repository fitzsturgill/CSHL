function setBoxMeanCoordinates(channel)
	global state
	global fsBoxWaveStruct
	
	
	
	figure(channel);

	k = waitforbuttonpress;

	if isempty(findobj(gcf, 'Type', 'axes'))
		disp('*** NO axes***');
		return
	end
		
		
	point1 = get(gca,'CurrentPoint');    % button down detected
	finalRect = rbbox;                   % return figure units
% 	set(gh.pcellControl.selectBoxButton, 'ForeGroundColor', [0 0 0]);
	setStatusString('');

	point2 = get(gca,'CurrentPoint');    % button up detected
	point1 = point1(1,1:2);              % extract x and y
	point2 = point2(1,1:2);
	p1 = min(point1,point2);             % calculate locations
	offset = abs(point1-point2);         % and dimensions
	x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
	y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
	x=round(x);
	y=round(y);
	fsBoxWaveStruct.x0=x(1);
	fsBoxWaveStruct.x1=x(2);
	fsBoxWaveStruct.y0=y(1);
	fsBoxWaveStruct.y1=y(3);
	fsBoxWaveStruct.channel=channel;
	
	waveO('fsBoxWave', [])
	
	