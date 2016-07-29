function out=openAndLoadUserSettings(fname)
% Allows user to select a settings file (*.ini) from disk and loads it
% Author: Bernardo Sabatini
	out=0;

	global state
	setStatusString('Loading user settings...');

	if nargin<1
		[fname, pname]=uigetfile('*.usr', 'Choose user settings file to load');
		if ~isnumeric(fname)
			periods=findstr(fname, '.');
			if any(periods)								
				fname=fname(1:periods(1)-1);
			else
				disp('openAndLoadUserSettings: Error: found file name without extension');
				setStatusString('Can''t open file...');
				return
			end		
			openusr(fullfile(pname, [fname '.usr']));
	
		else
			return
		end
	else
		openusr(fname);
	end	
	
	timerCallPackageFunctions('UserSettings');

	global gh	% BSMOD added 1/30/1 with lines below

	wins=fieldnames(gh);

	for winCount=1:length(wins)
		winName=wins{winCount};
		if isfield(state.windowPositions, [winName '_position']) 
			if length(getfield(state.windowPositions, [winName '_position']))==4
				oldPos=get(getfield(getfield(gh, winName), 'figure1'), 'Position');
				newPos=getfield(state.windowPositions, [winName '_position']);
				newPos(3)=oldPos(3);
				newPos(4)=oldPos(4);
				set(getfield(getfield(gh, winName), 'figure1'), 'Position', newPos);
			end
		end
	end
			
	setStatusString('');
