
function out=MP285FinishMove(check)
% MP285SetPos controls the position of the MP285
% 
% MP285SetPos 
% 
% Class Support
%   -------------
%   The input variable [x y z] contains the absolute motor target positions in microns. 
%   The optional paramter 'resolution' contains the resolution in nm (nanometers)
%	The value used depends on the MP285 microcode 
%		
%   Karel Svoboda 8/28/00 Matlab 6.0R
%	svoboda@cshl.org
% 	Modified 2/5/1 by Bernardo Sabatini to support global state and preset serialPortHandle
	out=0;
	global state
	if state.motor.motorOn==0
		return
	end

	if nargin<1
		check=1;
	end
	if ~state.motor.movePending
		disp('MP285FinishMove:  Error: Called with no move pending');
		out=1;
		return
	end
	
	status=state.internal.statusString;
	state.motor.movePending=0;
	if length(state.motor.serialPortHandle) == 0
		disp(['MP285SetPos: MP285 not configured']);
		state.motor.lastPositionRead=[];
		out=1;
		return;
	end

	try
		n=get(state.motor.serialPortHandle,'BytesAvailable');
		setStatusString('Waiting for move...');
		while n==0
			n=get(state.motor.serialPortHandle,'BytesAvailable');
		end
		temp=fread(state.motor.serialPortHandle,n); 
		if temp(1)~=13
			disp('MP285FinishMove: Error: CR not returned by MP285');
			out=1;
			return
		end
	catch
		disp('MP285FinishMove: Error in MP285 communication');
		out=1;
		return
	end		
			
	% check if position was attained
	if check
		setStatusString('Checking move...');
		xyzN=MP285GetPos;
		state.motor.lastPositionRead=xyzN;

		if fix(state.motor.requestedPosition*10) ~= fix(xyzN*10);
			setStatusString('Bad move.');
			disp(['MP285SetPos: Requested position not attained; check hardware']);
			state.motor.lastPositionRead=[];
			out=1;
			return
		end
	end
	
	state.motor.requestedPosition=[];
	state.motor.movePending=0;
	setStatusString(status);
	out=0;
