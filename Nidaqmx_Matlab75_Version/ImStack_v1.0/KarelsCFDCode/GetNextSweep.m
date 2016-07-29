function out=GetNextSweep(directory, base, sweep, extension, sign)
% GetNextSweep analyzes Directory\BaseSweep.Extension an returns Sweep+1
%		GetNextSweep returns a single number, the number following Sweep 
%		that satisfy the convention 'Base###.Extension'
% 		Note: extensions are of type '.ext'; Sign<0 means out<Sweep and vice versa
%		
%   Karel Svoboda 1/16/00 Matlab 5.3
%	 svoboda@cshl.org
%
%

if nargin < 4 
     disp(['-------------------------------']);  
     disp([' GetNextSweep v',version])
     disp(['-------------------------------']);
     disp([' usage: GetNextSweep(directory, base, sweep, extension, sign)']);
     disp([' usage: GetNextSweep(directory, base, sweep, extension)']);
     error(['### not enough parameters to proceed']); 
  end
  
  if nargin < 5 
     sign=1;
  end
  
% extension
%get list of existing sweepa 
sweeplist=GetSweepList(directory, base, extension);

i=1;

while (sweeplist(i) < sweep); i=i+1; end
'i'
sweeplist(i);
sweep;

%the following condition is required for the case where sweep is not
%a memeber of sweeplist
if (sign > 0) & (sweeplist(i)>sweep); i=i-1; end


if (sign < 0) & (i > 1); i=i-1; end
if (sign < 0) & (i == 1); 
   i=1; 
   disp(['lowest sweep reached']);
end
if (sign > 0) & (i < length(sweeplist)); i=i+1; end
if (sign > 0) & (i == length(sweeplist)) 
   disp(['largest sweep reached']);
end

out=sweeplist(i);