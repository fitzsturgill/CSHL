function out = percentile(a,f)
    % a = input matrix
    % returns value in a nearest to fractional percentile f (e.g. f = 0.9
    % for 90% percentile)
    a = a(~isnan(a)); % ignore NaNs
    if ~isempty(a)
        n = numel(a);
        a = reshape(a, n, 1);
        b = sort(a);
        out = b(min(n, max(1, round(n * f))));
    else
        out = NaN;
    end
    
    
%     nd = ndims(a);   
%     for d = nd:-1:1
%         a = sort(a, d);
%     end
    