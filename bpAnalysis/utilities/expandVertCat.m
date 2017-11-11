function data = expandVertCat(data, add, align)
    % to vertically concatenate data with expansion in case new data is
    % larger than old (in terms of number of columns)
    if nargin < 3
        align = 'left';
    end

    newColumns = max(0, size(add, 2) - size(data, 2));
    newRows = size(add, 1);
    switch align
        case 'left'
                data = [data NaN(size(data, 1), newColumns)];
                data = [data; NaN(newRows, size(data, 2))];
                data(end - newRows + 1:end, 1:size(add, 2)) = add;
        case 'right'
                data = [NaN(size(data, 1), newColumns) data];
                data = [data; NaN(newRows, size(data, 2))];
                data(end - newRows + 1:end, end - size(add, 2) + 1:end) = add;
    end
    