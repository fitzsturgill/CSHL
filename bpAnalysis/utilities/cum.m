function varargout = cum(a)
% function [sorted index]=cum(a) OR out = cum(a),  out being structure with
% fields sorted and index
    % for generation of cumulative probability plots
    
    if numel(a) ~= length(a)
        disp('error in cum: input argument not a row or column vector');
        return
    end
    
    a = a(~isnan(a));
    
    sorted = sort(a);
    index = (0:numel(a) - 1) / (numel(a) - 1);
    
    if nargout == 1
        out = struct('sorted', sorted, 'index', index);
        varargout = {out};
    elseif nargout == 2
        varargout = {sorted index};
    end
    