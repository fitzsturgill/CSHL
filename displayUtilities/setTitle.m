function setTitle(hAxes,string,fontsize,yoffset)
%
% INPUT
%   hAxes:
%   string:
%   fontsize:
%   position:
%

% Created: 5/17/10 - SRO

if nargin < 3
    fontsize = 7;
    yoffset = 0;
end

if nargin < 4
    yoffset = 0;
end

% Get handles
hTitle = get(hAxes,'Title');

% Get current position
pos = get(hTitle,'Position');
xpos = get(hAxes,'XLim');
pos(1) = mean(xpos);
ymax = max(get(hAxes,'YLim'));


pos(2) = ymax*(1 + yoffset);


% Set properties
set(hTitle,'String',string,'Fontsize',fontsize,'Position',pos,...
    'Color',[0 0 0], 'Interpreter', 'none');


