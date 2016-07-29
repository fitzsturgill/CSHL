function updatesAPLVariables(handle)
global state

% updatesAPLVariables.m******
% Callback function that updates all the parameters that need to be updatred whenever
% the samplesAcquiredPerLine variable changes.
%
% Written By: Thomas Pologruto
% Cold Spring Harbor Labs
% January 2, 2001

updatebinFactor;
