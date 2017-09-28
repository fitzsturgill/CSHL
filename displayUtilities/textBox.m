function textBox(s, h, position, fs)
    % make a text box in current axes, center and top justified
    % cell array of strings results in multiple lines
    if nargin < 2 || isempty(h)
        h = gca;
    end
    
    if nargin < 3 || isempty(position)
        position = [0.5 0.95];
    end
    
    if nargin < 4
        fs = 12;
    end
        
    
    text(position(1), position(2), s, ...
        'Units', 'normalized',...
        'Interpreter', 'none',...
        'HorizontalAlignment', 'center',...
        'FontSize', fs,...
        'VerticalAlignment', 'top'...
        );