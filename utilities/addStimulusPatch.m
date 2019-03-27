function varargout = addStimulusPatch(hAxes,position,str,color, alpha)
% function varargout = addStimulusBar(hAxes,position,str,color,linewidth)
% 
% INPUTS
%   hAxes(i): handles to axis 
%   position: [left right bottom top] left,right = x units,   bottom, top
%   (normalized 0 = bottom of axes, top = top
%   of y axes
%   str: String to be displayed above bar
%   color: 3-element RGB vector
%   linewidth: Thickness of line
%
% OUTPUTS
%   hLine: Handle to line object.

%   Created: 3/16/10 - SRO

if nargin < 3
    str = '';
    color = [0.5 0.5 0.5];
end

if nargin < 4
    color = [0.5 0.5 0.5];
end

if nargin < 5
    alpha = 0.5;
end

if length(position) == 2
    position(3) = 0;
end

if length(position) == 3
    position(4) = 1;
end



for i = 1:length(hAxes)

l = position(1);
r = position(2);
% y = position(3)-position(3)*0.03;
YLim = get(hAxes(i), 'YLim');
% b = Ylim(1) + (Ylim(2) - Ylim(1))/100;
% t = Ylim(2) - (Ylim(2) - Ylim(1))/100;
b = YLim(1) + diff(YLim) * max(position(3), 0.01);
t = YLim(1) + diff(YLim) * min(position(4), 0.99);
hLine = patch('XData', [l r r l], 'YData', [b b t t],'Parent',hAxes(i),'FaceColor',color, 'EdgeColor', 'none', 'FaceAlpha', alpha);
hAxes(i);
hText = text(mean([l r]),min(t + 0.1 * diff(YLim), YLim(1) + diff(YLim) * 0.9),str,'FontSize',8,'Color',color,...
    'HorizontalAlignment','center','Parent',hAxes(i));

end


% Outputs
varargout{1} = hLine;
% varargout{2} = hText;

