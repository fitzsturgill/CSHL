
function xyz=MP285GetPos
% MP285GetPos retrieves the position information from the MP285 controller
% 
% MP285GetPos 
% 
% Class Support
%   -------------
%   
%	the output [x y z] is the position of the MP285 in microns
%		
%   Karel Svoboda 8/28/00 Matlab 6.0R
%	 svoboda@cshl.org
% 	Modified 2/5/1 by Bernardo Sabatini to support global state and preset serialPortHandle

xyz=[];

global state
if state.motor.motorOn==0
	return
end
	
if length(state.motor.serialPortHandle) == 0
	disp(['MP285GetPos: MP285 not configured.']);
	xyz=[];
	state.motor.lastPositionRead=[];
	return
end

%whos state.motor.serialPortHandle;
% get all the junk out
MP285Flush;

mp285Error=0;

% read position
try
	fwrite(state.motor.serialPortHandle, [99 13]); 		%'c'CR
catch
	mp285Error=1;
end

if ~mp285Error
	try
		array = fread(state.motor.serialPortHandle, 3, 'long');		% read position information (12bytes) including CR (1 byte)
	catch
		mp285Error=1;
	end
end

if ~mp285Error
	try
		dummy = fread(state.motor.serialPortHandle, 1);		% read position information (12bytes) including CR (1 byte)
	catch
		mp285Error=1;
	end
end

if mp285Error
	disp('mp285GetPos: Error in MP285 Communication');
	disp(lasterr)
	setStatusString('MP285 Error. Reset?');
	state.motor.lastPositionRead=[];
	return
end


if length(array)<3 | length(dummy)<1				% check if data is avaiable
	disp(['MP285GetPos: MP285 position data not available ']);
	xyz=[];
	state.motor.lastPositionRead=[];
	return;	
end
xyz=reshape(array,1,3)./[state.motor.calibrationFactorX state.motor.calibrationFactorY state.motor.calibrationFactorZ]/10;

state.motor.lastPositionRead=xyz;

