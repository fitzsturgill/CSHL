function varargout = addStimulusPatch(hAxes,position,str,color)
% function varargout = addStimulusBar(hAxes,position,str,color,linewidth)
% 
% INPUTS
%   hAxes(i): handles to axis 
%   position: [left right]
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



for i = 1:length(hAxes)

l = position(1);
r = position(2);
% y = position(3)-position(3)*0.03;]
Ylim = get(hAxes(i), 'YLim');
b = Ylim(1) + (Ylim(2) - Ylim(1))/100;
t = Ylim(2) - (Ylim(2) - Ylim(1))/100;
hLine = patch('XData', [l r r l], 'YData', [b b t t],'Parent',hAxes(i),'FaceColor',color, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hAxes(i);
% hText = text(mean([l r]),y+0.08*y,str,'FontSize',8,'Color',color,...
%     'HorizontalAlignment','center','Parent',hAxes(i));

end


% Outputs
varargout{1} = hLine;
% varargout{2} = hText;

