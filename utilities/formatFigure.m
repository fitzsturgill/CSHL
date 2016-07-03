function formatFigure(aspect)
    if nargin < 1
        aspect = [2 1];        
    end
    scaleFactor = 1;
    fontSize = 10;
    printWidth = 10; 
    paperUnits = 'centimeters';
    screenPosition = [100 100 aspect * 100 * scaleFactor];
    paperPosition = [0 0 aspect * 5 * scaleFactor];
    % in points 8.5 x 11 is 595 x 770
    
    set(gca, 'TickDir', 'out');
    set(gcf, 'PaperUnits', paperUnits);
    set(gcf, 'Position', screenPosition, 'PaperPosition', paperPosition);
    