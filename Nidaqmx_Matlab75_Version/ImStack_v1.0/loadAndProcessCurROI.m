function loadAndProcessCurROI
% This will laod the current saved ROI and also process it according to
% computeTimeSeriesMean protocol.

flag=loadCurrentROI;
if flag
    computeTimeSeriesMean;
end