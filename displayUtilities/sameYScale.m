function sameYScale(ha)
% for a series of axes, set y axes limits to the lowest and highest common
% values


    ylims = get(ha, 'YLim');
    ymin = min(cellfun(@(x) x(1), ylims)); 
    ymax = max(cellfun(@(x) x(2), ylims));
    
    set(ha, 'YLim', [ymin ymax]);