function [sorted index]=cum(a)
    % for generation of cumulative probability plots
    
    if numel(a) ~= length(a)
        disp('error in cum: input argument not a row or column vector');
        return
    end
    
    a = a(~isnan(a));
    
    sorted = sort(a);
    index = (0:numel(a) - 1) / (numel(a) - 1);
    