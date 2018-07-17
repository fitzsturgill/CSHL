function inds = sort_inds(A,dim,mode)
% Returns the indices of a sorting of x.
%
% Equivalent to 
%       [~,inds] = sort(x,ascend_descend)
%       return inds
%
% (C) Alex Vaughan, March 2015q

assert(isnumeric(A),'Inputs must be numeric')

if nargin < 3,
    mode = 'ascend';
end

% By default, sort on the first dimension whose size > 1.
if nargin < 2,
    dim = 1;
    for this_dim = 1:ndims(A),
       if size(A,this_dim) > 1,
          dim = this_dim;
       end
    end
end

[~,inds] = sort(A,dim,mode);

