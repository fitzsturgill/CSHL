function zoomValues(axis)

XLimits = round(get(axis, 'XLim'))
YLimits = round(get(axis, 'YLim'))

% Linear relationship beween pixels in each directiona and 
% intensities in each direction.

slopeX = 
intensity