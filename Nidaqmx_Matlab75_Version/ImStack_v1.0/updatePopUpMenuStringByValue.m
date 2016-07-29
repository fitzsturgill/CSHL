function updatePopUpMenuStringByValue(object, value, string)
global state gh

% this function will update a popupmenu string by the value entered
% at the position desired.

a = get(object, 'String');
if iscell(a)
	size = size(a,1);
else
	size=1;
end

if value > size
	disp('Error: value larger than string. Use GUIDE to increase String Size.');
else
	if ischar(a)		
		set(object, 'String', string);
	else
		a{value} = string;
		set(object, 'String', string);
	end
end
