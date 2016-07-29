function zoomIn(axis)
global gh state

rectArea = getrect(axis);
Xlim = [rectArea(1) (rectArea(1)+rectArea(3))];
Ylim = [rectArea(2) (rectArea(4)+rectArea(2))];
set(axis, 'Xlim', Xlim, 'Ylim', Ylim);
