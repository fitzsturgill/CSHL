function output = to_cell(input)
% Converts ~any input to a cell.
% Alex Vaughan, 2015

    input = input(:);
    if strcmp(class(input),'cell')
        output = input;
    elseif isnumeric(input)
        output = num2cell(input(:));
    elseif ischar(input)
        output = {input};
    elseif isstruct(input)
        for i = 1:length(input)
            output{i} = input(i);
        end
    end
end