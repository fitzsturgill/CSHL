function updateNextPopUpWithString(object, string)
global state gh

% this function will update a popupmenu string in the next value 
	b = get(object, 'String');
	
	if ischar(b)
		a{1} = b;	
		sizeCell = length(a);
		if sizeCell == 2 & strcmp(a{2}, '')	
			newvalue = sizeCell;
			a{newvalue} = string;
			set(object, 'String', a);
		else
			newvalue = sizeCell + 1;
			a{newvalue} = string;
			set(object, 'String', a);
		end
	else
		a = b;
		sizeCell = length(b);
		if sizeCell == 2 & strcmp(a{2}, '')	
			newvalue = sizeCell;
			a{newvalue} = string;
			set(object, 'String', a);
		else
			newvalue = sizeCell+1;
			a{newvalue} = string;
			set(object, 'String', a);
		end
	end
	
			