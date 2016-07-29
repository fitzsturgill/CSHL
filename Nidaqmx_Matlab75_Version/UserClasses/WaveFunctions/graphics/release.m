function release(axishandle,varargin)
% Function releases scale and the resacaleAxis function will operate on it.
% USe freeze(ax) to forbade autoscale 
%
% Sets the axis UserData.autoscale=1 on the ax.
%
% Operates on current axes by default.
% 
% varargin in can now be 'autoscale' which can be set to 2 for Matlab
% scaling...

if nargin == 0	% No inputs...look for GCA
	a=findobj('Type', 'axes');  % Are there any axes?
	if ~isempty(a)              
		ax=gca;
	else
		error('release: No axes defined.');
	end
elseif nargin>=1	% axis handle...check
	if istype(axishandle,'axes')
		ax=axishandle;
	else
		error('release: Bad Axes handle.');
	end
end

autoscale=1;
% Parse input parameter pairs and rewrite values.
counter=1;
while counter+1 <= length(varargin)
	eval([varargin{counter} '=' num2str(varargin{counter+1}) ';']);
	counter=counter+2;
end

userData = get(ax, 'UserData'); % Check if frozen....the dont autoscale.
if isstruct(userData)	% Check format...
	if isfield(userData, 'autoScale')
		userData=setfield(userData, 'autoScale', autoscale); 	%Dont AutoScale
		set(ax,  'UserData', userData);
		rescaleAxis(ax);
	else
		return	
	end
else
	error('release: axis UserData in incorrect format.');
end


