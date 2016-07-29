function setXYsameLimit(hax, includeOrgin)
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

all = [xlim ylim];
minL = min(all);
maxL = max(all);
if includeOrgin
    lim = [0 maxL];
else
    lim = [minL maxL];
end

set(hax,'Xlim',lim,'YLim',lim);