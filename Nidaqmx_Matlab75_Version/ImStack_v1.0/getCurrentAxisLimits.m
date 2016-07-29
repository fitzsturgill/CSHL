function [XStart, XEnd, YStart, YEnd] = getCurrentAxisLimits(axis)


Xlimit = get(axis, 'XLim');
Ylimit = get(axis, 'YLim');

XStart = ceil(Xlimit(1,1));
if XStart < 1
	XStart =1;
end

XEnd = fix(Xlimit(1,2));
YStart = ceil(Ylimit(1,1));
if YStart < 1
	YStart =1;
end

YEnd = fix(Ylimit(1,2));

