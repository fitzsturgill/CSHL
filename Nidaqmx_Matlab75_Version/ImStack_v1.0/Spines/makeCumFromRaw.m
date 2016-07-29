function out = makeCumFromRaw(input)
% Make a column of a histogram to a normalized cumulative distribution

if ndims(input) == 1
	
	Sum = sum(input);
	input = input./Sum;
	for i = 1:size(input,1)
		out(i) = input(i)+ sum(input(1:i-1));
	end

elseif ndims(input) == 2
	
	for j = 1:size(input,2)
		Sum = sum(input(:,j));
		input(:,j) = input(:,j)./Sum;
		for i = 1:size(input(:,j),1)
			out(i,j) = input(i,j)+ sum(input(1:i-1,j));
		end
	end
else
	beep;
	disp('Dimension of input must be <= 2');
	return
end

		
