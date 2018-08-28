function cell_grid = ndgrid_cell(varargin)
% Generate ndgrid-style combinations of inputs, modified to allow for any
% input type (array, cell array, struct) while returning cell arrays.
%
% Inputs:
%   varargin  ::    ~any n-d array object to be iterated over.
%                   Iterates over cell arrays, numeric arrays, and
%                   n-dimensional structs.  Char arrays are wrapped as a cell
%                   array with one element.
%
% Outputs:
%   cell_grid ::    a cell array of size nargin, with each element
%                   representing a cell array ndgrid across the elements 
%                   of an inputThis is roughtly similar to indexing the 
%                   varargin inputs with the output of an indexed ndgrid.
%
%%
% Notes:
%   This function builds on the answer here: 
%        http://stackoverflow.com/questions/8492277/matlab-combinations-of-an-arbitrary-number-of-cell-arrays
%   The original code read 
%       indices = fliplr(indices), 
%       indices{:} = ndgrid(indices{:})
%       cellfun(...,fliplr(indices))
%   this ensures that the first input is the most significant bit (ie., changes most slowly).
%   This behavior has been changed to a '--sorted' flag.
%
% Alex Vaughan, 2015

do_sort = false;
if ischar(varargin{end}) && strcmp(varargin{end},'--sorted')
   do_sort = true;
   varargin = fliplr({varargin{1:end-1}});
end

sizeVec = cellfun('prodofsize', varargin);
indices = arrayfun(@(n) {1:n}, sizeVec);
[indices{:}] = ndgrid(indices{:});

% Ensure that all inputs are cell arrays
for i = 1:length(varargin),
    varargin{i} = to_cell(varargin{i});
end

% Original code from stackoverflow.
% Returns a single cell array, with one column per input variable and one
% row per combination of variables.
%combMat = cellfun(@(c,i) {reshape(c(i(:)), [], 1)}, varargin, fliplr(indices));
%combMat = [combMat{:}];

% Return cell arrays with the original size of ndgrid.
cell_grid = cellfun(@(c,i) {c(i)}, varargin, (indices));

if do_sort,
   cell_grid = fliplr(cell_grid);
end

end


