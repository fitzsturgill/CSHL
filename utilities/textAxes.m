function [ax txt] = textAxes(h, txt, fsize)
% create an empty axes to display only text on a figure
% h - figure handle
% text (optional) text to display
    
    if nargin < 3
        fsize =  8;
    end
    
    if isa(h, 'matlab.ui.Figure')
        ax = axes('Parent', h);
        axis off
    else
        axes(h);
        cla;
        axis off;
    end
    

        
    
    if nargin > 1
        x=0;
        y=0;
        txt = text(...
            'String', txt,...
            'Units', 'normalized',...
            'Position', [.5 .5],...
            'FontSize', fsize,...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'middle',...
            'Interpreter', 'none'...
            );
    else
        txt = [];
    end