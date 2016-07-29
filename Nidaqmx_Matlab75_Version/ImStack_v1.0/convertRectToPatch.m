function patchData=convertRectToPatch(rectInfo)
%This function takes the results of a call to getrect
% anfd formats it to look like a patch object polygon

if nargin ~= 1
    error('convertRectToPatch: Only 1 1x4 array allowed as input');
end

if length(rectInfo) ~= 4
     error('convertRectToPatch: Only 1 1x4 array allowed as input');
end

patchData=[rectInfo(1) rectInfo(1) rectInfo(1)+rectInfo(3) rectInfo(1)+rectInfo(3); ...
        rectInfo(2) rectInfo(2)+rectInfo(4) rectInfo(2)+rectInfo(4) rectInfo(2)];



    
