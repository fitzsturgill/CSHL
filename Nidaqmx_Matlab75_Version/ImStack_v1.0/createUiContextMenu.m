function cmenuHandle=createUiContextMenu(varargin)
% This function creates a uicontextMenu and passes the handle as an output.
% The inputs are strings containing the fucntions to be executed
% They are param-function string pairs, where the param is what will be displayed 
% in the menu, and the function is the fucntion to be executed

if nargin < 2
    error('createUiContextMenu: Must pass in name-function string pairs as inputs');
end
cmenuHandle = uicontextmenu;
paramsInput=varargin;
while length(paramsInput) >= 2
    uimenu(cmenuHandle, 'Label', paramsInput{1}, 'Callback', paramsInput{2});
    paramsInput=paramsInput(3:end);
end
