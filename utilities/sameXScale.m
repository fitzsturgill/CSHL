function sameXScale(ha)
% for a series of axes, set x axes limits to the lowest and highest common
% values


    xlims = get(ha, 'XLim');
    xmin = min(cellfun(@(x) x(1), xlims)); 
    xmax = max(cellfun(@(x) x(2), xlims));
    
    set(ha, 'XLim', [xmin xmax]);