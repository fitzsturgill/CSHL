function formatFigurePoster(wh, h)
%% note for .pdf paper units matter (or other paged formats)
%% for .eps, .tif, .jpg,  screen units matter
    if nargin < 1
        wh = [4 3]; % to match [4 3] aspect ratio
    end
    
    if nargin < 2
        h = gcf;
    end
    aspect = wh(2) / wh(1);  % normalize by width (make width = 1)
    units = 'inches';
    screenPosition = [1 1 wh(1) wh(2)];
    paperPosition = [0 0 wh(1) wh(2)];
    % in points 8.5 x 11 is 595 x 770
    
    ax = findobj(h, 'Type', 'Axes');
    set(ax, 'TickDir', 'out', 'LineWidth', 1, 'FontSize', 24, 'FontName', 'Arial', 'Box', 'off');
    set(gcf, 'Units', units, 'PaperUnits', units, 'PaperSize', [wh(1) wh(2)],...
        'ToolBar', 'none', 'MenuBar', 'none',  'DockControls', 'off');
    set(gcf, 'Position', screenPosition, 'PaperPosition', paperPosition);