function sameXYScale(ha)
% for a series of axes, set x and y axes limits to the lowest and highest common
% values


    ylims = get(ha, 'YLim');
    ymin = min(cellfun(@(x) x(1), ylims)); 
    ymax = max(cellfun(@(x) x(2), ylims));
    xlims = get(ha, 'XLim');
    xmin = min(cellfun(@(x) x(1), xlims)); 
    xmax = max(cellfun(@(x) x(2), xlims));
    
    xymin = min(xmin, ymin);
    xymax = max(xmax, ymax);
    set(ha, 'XLim', [xymin xymax], 'YLim', [xymin xymax]);