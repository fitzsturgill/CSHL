function A_reordered = reorder_by_sort(A,B,dim,mode)
% Returns the entries of A reordered by the sort index of B
%
% Equivalent to 
%       [~,inds] = sort(B,ascend_descend)
%       A_reordered = inds(A)
%       return A_reordered
%
% (C) Alex Vaughan, March 2015

assert isnumeric(x)
assert ndims(squeeze(A))==1

if nargin < 3,
    mode = 'ascend';
end

% By default, sort on the first dimension whose size > 1.
if nargin < 3,
    dim = 1;
    for this_dim = 1:ndims(B),
       if size(B,this_dim) > 1,
          dim = this_dim;
       end
    end
end

[~,reorder_inds] = sort(B,dim,mode);
A_reordered = reorder_inds(A)

