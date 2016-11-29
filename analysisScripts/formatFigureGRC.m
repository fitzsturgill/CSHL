function formatFigureGRC(wh)
%% note for .pdf paper units matter (or other paged formats)
%% for .eps, .tif, .jpg,  screen units matter
    if nargin < 1
        wh = [93 42];
    end
    aspect = wh(2) / wh(1);  % normalize by width (make width = 1)
    units = 'centimeters';
    screenPosition = [5 5 wh(1) wh(2)];
    paperPosition = [0 0 wh(1) wh(2)];
    % in points 8.5 x 11 is 595 x 770
    
    set(gca, 'TickDir', 'out', 'LineWidth', 1);
    set(gcf, 'Units', units, 'PaperUnits', units);
    set(gcf, 'Position', screenPosition, 'PaperPosition', paperPosition);
    end