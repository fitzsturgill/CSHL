
function out=MP285SetPos(xyz, resolution, checkPosition)
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

out=1;
global state
if state.motor.motorOn==0
	return
end

if nargin < 1
     disp(['-------------------------------']);  
     disp([' MP285SetPos v',version])
     disp(['-------------------------------']);
     disp([' usage: MMP285SetPos([x y z])']);
     error(['### incomplete parameters; cannot proceed']); 
end 

if nargin < 2
     resolution=state.motor.resolution;
end

if nargin < 3
	checkPosition=1;
end

if isempty(resolution)
     resolution=state.motor.resolution;
end

if isempty(checkPosition)
	checkPosition=1;
end

if length(xyz) ~=3
     disp(['-------------------------------']);  
     disp([' MP285SetPos v',version])
     disp(['-------------------------------']);
     disp([' usage: MP285SetPos([x y z])'])
     error(['### incomplete or ambiguous parameters; cannot proceed']); 
end 

if length(state.motor.serialPortHandle) == 0
	disp(['MP285SetPos: MP285 not configured']);
	state.motor.lastPositionRead=[];
	out=1;
	return;
end
 
% convert microns to units of nm  mod resolution (i.e. 100nm resolution);
xyz2=fix(xyz*state.motor.resolution).*	...
	[state.motor.calibrationFactorX state.motor.calibrationFactorY state.motor.calibrationFactorZ];

% flush all the junk out
MP285Flush;

% temp=MP285Comp14ByteArr(xyz);
try
	fwrite(state.motor.serialPortHandle, 'm');
	fwrite(state.motor.serialPortHandle, xyz2, 'long');
	fwrite(state.motor.serialPortHandle, 13);
	out=fread(state.motor.serialPortHandle,1);
catch
	disp(['MP285SetPos: MP285 communication eror.']);
	return
end

if out ~= 13; 
	disp(['MP285SetPos: MP285 return an error.  Unsure of movement status.']); 
	MP285Flush;
	state.motor.lastPositionRead=[];
	out=1;
	return;
end				% check if CR was returned

% check if position was attained
if checkPosition
	xyzN=MP285GetPos;
	if isempty(xyzN)
		disp(['MP285SetPos: Unable to check movement.']);
	elseif fix(xyz*state.motor.resolution) ~= fix(xyzN*state.motor.resolution); 
		disp(['MP285SetPos: Requested position not attained; check hardware']);
		state.motor.lastPositionRead=[];
		out=1;
		return;
	end		
end

out=0;



