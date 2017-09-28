function removeAxesLabels(hAxes,b)
% function removeAxesLabels(hAxes,b)
% INPUT:
%   hAxes: Handles to axis or vector of axes b: Struct with flags for
%   removing specific labels .rmXlabel, .rmYlabel, .rmXticklabel,
%   .rmYticklabel

% Created: 3/10 - SRO
% Modified: 10/20/10 - SRO

if nargin < 2
    b.xl = 1;
    b.yl = 1;
    b.xtl = 1;
    b.ytl = 1;
end

for i = 1:length(hAxes)
    if b.xl
        xlabel(hAxes(i),'');
    end
    
    if b.yl
         ylabel(hAxes(i),'');
    end
    
    if b.xtl
        set(hAxes(i),'XTickLabel',{});
    end
 
    if b.ytl
        set(hAxes(i),'YTickLabel',{});
    end                  
end

