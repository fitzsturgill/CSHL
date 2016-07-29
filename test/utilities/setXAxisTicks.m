function setXAxisTicks(hAxis)
% 
% Adds 2 ticks to x-axis at 25% and 75% of yLim.
%
%

% Created: 8/13/10 - SRO

tempXLim = get(hAxis,'XLim');
adj = 0.25*diff(tempXLim);
tempXLim = tempXLim + [adj -adj];


if max(abs(tempXLim)) > 10
    tempXLim = roundn(tempXLim,0);
elseif max(abs(tempXLim)) > 1
    tempXLim = roundn(tempXLim,-1);
elseif max(abs(tempXLim)) > 0.1
    tempXLim = roundn(tempXLim,-2);
else
    tempXLim = roundn(tempXLim,-3);
end

if diff(tempXLim)
    set(hAxis,'XTick',tempXLim,'XTickLabel',num2cell(tempXLim));
elseif ~diff(tempXLim)
    tempXLim = tempXLim + [-1.1 1.1].*tempXLim;
    if any(tempXLim)
        set(hAxis,'XTick',tempXLim,'XTickLabel',num2cell(tempXLim))
    else
        set(hAxis,'XTick',tempXLim(1),'XTickLabel',num2cell(tempXLim(1)))
    end
end