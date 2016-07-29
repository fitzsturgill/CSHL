function updatefillFractionVariables(handle)
global state

% updatefillFractionVariables.m******
% Callback function that updates all the parameters that need to be updatred whenever
% the fillFraction variable changes.
%
% Written By: Thomas Pologruto
% Cold Spring Harbor Labs
% January 14, 2001

updatePixelTime;
setAcquisitionParameters;
updatebinFactor;