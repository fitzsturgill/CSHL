function cellOutput = removeFieldFromCellArray(cellInput, indextoremove, row)
global state gh

% removes the field given by indextoremove from cell arrays 
% row si 1 if the cell array is a row array, and it is 0 if it is a column cell array

numberofcellsfinal = 1;
lengthInput = length(cellInput);

if lengthInput > 1
	
	fieldcounter = lengthInput-1;
	if row
		for counter = 1:fieldcounter
			if counter == indextoremove
				numberofcellsfinal = numberofcellsfinal +1;
				cellOutput{counter} = cellInput{numberofcellsfinal};
			else 
				cellOutput{counter} = cellInput{numberofcellsfinal};	
			end
			numberofcellsfinal = numberofcellsfinal +1;
		end

	else
		for counter = 1:fieldcounter
			if counter == indextoremove
				numberofcellsfinal = numberofcellsfinal +1;
				cellOutput{counter,1} = cellInput{numberofcellsfinal};
			else 
				cellOutput{counter,1} = cellInput{numberofcellsfinal};	
			end
			numberofcellsfinal = numberofcellsfinal +1;
		end
	end

else
	cellOutput = {};
end


			
			