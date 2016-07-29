function h = addUnityLine(hax,c)
% function h = addUnityLine(hax,c)
%
% INPUT
%   hax: Axes handle
%   c: 3-element rgb vector for color of unity line

% Created: SRO - 6/10/11

if nargin < 1 || isempty(hax)
    hax = gca;
end

if nargin < 2 || isempty(c)
    c = [0.7 0.7 0.7];
end


% Get axes limits
xlim = get(hax,'XLim');
ylim = get(hax,'YLim');

p1 = min([min(xlim) min(ylim)]);
p2 = min([max(xlim) max(ylim)]);

h = line('Parent',hax,'XData',[p1 p2],'YData',[p1 p2]);

set(h,'Color',c);

