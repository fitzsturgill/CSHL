function updatemsPerLineVariables(handle)
global state gh

% updatemsPerLineVariables.m******
% Callback function that updates all the parameters that need to be updatred whenever
% the msPerLine variable changes.
%
% Written By: Thomas Pologruto
% Cold Spring Harbor Labs
% January 2, 2001

updatePixelTime;
setAcquisitionParameters;
updatebinFactor;