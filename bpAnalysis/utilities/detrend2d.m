function out = detrend2d(data, dim)
    if nargin < 2
        dim = 1;
    end
    
    if dim == 2
        data = data';
    end
    out = zeros(size(data));
    
    for counter = 1:size(data, 1)
        out(counter, :) = detrend(data(counter, :));
    end
    
    if dim == 2
        out = out';
    end
        
    