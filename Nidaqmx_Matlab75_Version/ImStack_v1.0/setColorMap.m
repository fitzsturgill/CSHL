function setColorMap(map)
global gh state

switch map
case 'green'
	colormap(makeColorMap('green',8));
case 'red'
	colormap(makeColorMap('red',8));
case 'blue'
	colormap(makeColorMap('blue',8));
case 'gray'
	colormap(makeColorMap('gray', 8));
otherwise
	eval(['colormap(' map ');']);
end
