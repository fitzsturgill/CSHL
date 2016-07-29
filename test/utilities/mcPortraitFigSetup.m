function h = mcPortraitFigSetup(h)
%function h = mcPortraitFigSetup(h)
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
width = 8.5-leftMargin-rightMargin;
height = 11-topMargin-bottomMargin;
PaperPosition = [leftMargin bottomMargin width height];

% Set monitor position of figure
ppi = get(0, 'ScreenPixelsPerInch');
screenWidth = width * ppi;
screenHeight = height * ppi;
monitorPosition = [0 0 screenWidth screenHeight];


set(h,'PaperOrientation','portrait','PaperPosition',PaperPosition, ...
    'Position',monitorPosition)

