function superclear
% A stupid function that clears everything.
%
try,
evalin('caller','dbquit')
end
evalin('caller','clear')
evalin('caller','close all')
evalin('caller','clc')