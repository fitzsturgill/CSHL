function roiWindowBUpFcn
global gh state
set(gcf,'WindowButtonMotionFcn','');

% Do analysis
roiStats;