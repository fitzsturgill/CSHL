function h = mcLandscapeFigSetup(h)
% hFig = landscapeFigSetup(varargin)
%
%
%

% Created: 3/15/10 - SRO



if nargin==1
    h = figure(h);
else
    h = figure;
end


% Set paper position in inches (where figure will be printed on paper)
leftMargin = 0.15;
rightMargin = 0.15;
bottomMargin = .3;
topMargin = 0.4;
width = 11-leftMargin-rightMargin;
height = 8.5-topMargin-bottomMargin;
PaperPosition = [leftMargin bottomMargin width height];

% Set monitor position of figure
ppi = get(0, 'ScreenPixelsPerInch');
screenWidth = width * ppi;
screenHeight = height * ppi;
monitorPosition = [0 0 screenWidth screenHeight];


set(h,'PaperOrientation','landscape','PaperPosition',PaperPosition, ...
    'Position',monitorPosition)

