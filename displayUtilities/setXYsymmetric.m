function setXYsymmetric(hax, includeOrgin)
%
%
%
%

% Created: SRO - 6/16/11

if nargin < 1
    hax = gca;
end

if nargin < 2
    includeOrgin = 1;
end

xlim = get(hax,'XLim');
ylim = get(hax,'YLim');


maxX = max(abs(xlim));
maxY = max(abs(ylim));

set(hax,'Xlim',[-maxX maxX],'YLim',[-maxY maxY]);