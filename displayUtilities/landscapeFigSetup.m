function hFig = landscapeFigSetup(varargin)
% hFig = landscapeFigSetup(varargin)
%
%
%

% Created: 3/15/10 - SRO

if nargin > 0 
    hFig = varargin{1};
else
    hFig = [];
end

if ~isempty(hFig)
    hFig = figure(hFig,'Visible','off');
else
    hFig = figure('Visible','off');
end

% Set monitor position of figure
monitorPosition = [154 59 1123 742];

% Set paper position in inches (where figure will be printed on paper)
leftMargin = 0.15;
rightMargin = 0.15;
bottomMargin = 0;
topMargin = 0.4;
PaperPosition = [leftMargin bottomMargin (11-leftMargin-rightMargin) (8.5-topMargin-bottomMargin)];

set(hFig,'PaperOrientation','landscape','PaperPosition',PaperPosition, ...
    'Position',monitorPosition,'Color',[0.94 0.94 0.94])

% Get screen size
ScreenSize = get(0,'ScreenSize');

% Make figure visible
set(hFig,'Visible','on')


