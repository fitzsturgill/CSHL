function setLineColors(handle, color)

% Set lines in axis to color for better contrast
type = get(handle, 'Type');
if ~strcmp(type, 'axes')
	return
end

a = get(handle, 'Children');
switch color	
case 'black'
	col = [ 0 0 0];
case 'gray'
	col = [ .8 .8 .8];
case 'red'
	col = [ 1 0 0];
case 'green'
	col = [ 0 1 0];
case 'blue'
	col = [ 0 0 1];
otherwise
	col = [ 0 0 1];
end

for i = 1:length(a)
	lineT = get(a(i), 'Type');
	if strcmp(lineT,'line')
		set(a(i), 'Color', col);
	end
end
parent = get(handle, 'Parent');
set(parent, 'Color', [1 1 1]);
 