function out = getMenuIndex2(handle, name);
% get value associated with the input name

string = get(handle, 'String');
if ischar(string)
	if strcmp(string, name)
		out = 1;
	else 
		disp('Not Found');
		return
	end
elseif iscell(string)
	for i = 1:length(string)
		if strcmp(string(i), name)
			out = i;
		end
	end
end

	
	