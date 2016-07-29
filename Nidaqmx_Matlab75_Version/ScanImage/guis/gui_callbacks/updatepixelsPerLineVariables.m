function updatepixelsPerLineVariables(handle)
global state

% updatepixelsPerLineVariables.m******
% Callback function that updates all the parameters that need to be updatred whenever
% the pixelsPerLine variable changes.
%
% Written By: Thomas Pologruto
% Cold Spring Harbor Labs
% January 2, 2001

updatebinFactor;
updatePixelTime;
