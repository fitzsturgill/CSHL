function out = percentile(a,f, dim)
    % a = input matrix
    % returns value in a nearest to fractional percentile f (e.g. f = 0.9
    % for 90% percentile)
    if nargin < 3
        dim = 0;
    end
    if dim > 2
        error('can not handle more than 2 dimensions currently');
    end
    if ~dim
        out = getP(a, f);
    else
        
        if dim == 2
            out = NaN(size(a, 1), 1);
            a = a';
        else
            out = NaN(1, size(a, 2));
        end
        for counter = 1:size(a, 2)
            out(counter) = getP(a(:,counter), f);
        end
    end
end
        
    
        
function p = getP(a, f)
    a = a(~isnan(a)); % ignore NaNs
    if ~isempty(a)
        b = sort(a);
        n = length(a);
        p = b(min(n, max(1, round(n * f))));
    else
        p = NaN;
    end
end
    
    
%     nd = ndims(a);   
%     for d = nd:-1:1
%         a = sort(a, d);
%     end
    