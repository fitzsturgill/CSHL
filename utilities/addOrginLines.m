function h = addOrginLines(hax,c)
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

xp1 = [min(xlim) 0];
xp2 = [max(xlim) 0];

yp1 = [0 min(ylim)];
yp2 = [0 max(ylim)];

xh = line('Parent',hax,'XData',[xp1(1) xp2(1)],'YData',[xp1(2) xp2(2)]);
yh = line('Parent',hax,'XData',[yp1(1) yp2(1)],'YData',[yp1(2) yp2(2)]);

h = [xh yh];
set(h,'Color',c);
set(hax, 'XLim', xlim, 'YLim', ylim); % FS MOD 10/2018

