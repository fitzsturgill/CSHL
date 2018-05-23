function formatFigure(varargin)
    defaults = {...
        'aspect', [2 1];... % aspect ratio w x h (i believe)
        'scaleFactor', 1;...
        };

    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    
    
    fontSize = 10;
    printWidth = 10; 
    paperUnits = 'centimeters';
    screenPosition = [100 100 s.aspect * 100 * s.scaleFactor];
    paperPosition = [0 0 s.aspect * 5 * s.scaleFactor];
    % in points 8.5 x 11 is 595 x 770
    
    set(gca, 'TickDir', 'out');
    set(gcf, 'PaperUnits', paperUnits);
    set(gcf, 'Position', screenPosition, 'PaperPosition', paperPosition);
    