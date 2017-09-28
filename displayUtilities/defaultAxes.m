function defaultAxes(hAxes,xdta,ydta,fontsize)
% function defaultAxes(hAxes,xdta,ydta,fontsize)
%
% INPUTS
%   hAxes:
%   xdta:
%   ydata:
%   fontsize:

% Created: 3/15/10 - SRO
% Modified: 5/16/10 - SRO


if nargin < 2 || isempty(xdta)
    xdta = 0.10;    % Distance to axis
end

if nargin < 3 || isempty(ydta)
    ydta = 0.13;
end

if nargin < 4 || isempty(fontsize)
    fontsize = 8;
end

% Set color of axes
% axesGray(hAxes)

% Set tick out
set(hAxes,'TickDir','out');

% Set font size
set(hAxes,'FontSize',fontsize);

% Remove box
% box off

% Set position of tick label

% Set position and font size of axis label
% for i = 1:length(hAxes)
pta = 0.5;      % Parallel to axis
if length(hAxes) > 1
    set(cell2mat(get(hAxes,'XLabel')),'Units','normalized','Position',[pta -xdta 1],'FontSize',fontsize);
    set(cell2mat(get(hAxes,'YLabel')),'Units','normalized','Position',[-ydta pta 1],'FontSize',fontsize);
else
    set(get(hAxes,'XLabel'),'Units','normalized','Position',[pta -xdta 1],'FontSize',fontsize);
    set(get(hAxes,'YLabel'),'Units','normalized','Position',[-ydta pta 1],'FontSize',fontsize);
end

