function setAxisTicks(hAxis)
% 
% Adds 2 ticks to y-axis at 25% and 75% of yLim.
%
%

% Created: 8/13/10 - SRO

tempYLim = get(hAxis,'YLim');
adj = 0.25*diff(tempYLim);
tempYLim = tempYLim + [adj -adj];


if max(abs(tempYLim)) > 10
    tempYLim = roundn(tempYLim,0);
elseif max(abs(tempYLim)) > 1
    tempYLim = roundn(tempYLim,-1);
elseif max(abs(tempYLim)) > 0.1
    tempYLim = roundn(tempYLim,-2);
else
    tempYLim = roundn(tempYLim,-3);
end

if diff(tempYLim)
    set(hAxis,'YTick',tempYLim,'YTickLabel',num2cell(tempYLim));
elseif ~diff(tempYLim)
    tempYLim = tempYLim + [-1.1 1.1].*tempYLim;
    if any(tempYLim)
        set(hAxis,'YTick',tempYLim,'YTickLabel',num2cell(tempYLim))
    else
        set(hAxis,'YTick',tempYLim(1),'YTickLabel',num2cell(tempYLim(1)))
    end
end
