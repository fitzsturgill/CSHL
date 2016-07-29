function updateGUIByGlobalCell(fieldname, location)
global state gh

% this function will update the GUI's associated with a global variable when the 
% global is a cell array.

[structureName,structureToEnd,lastField] = structNameParts(fieldname);
eval([fieldname ' = ' structureToEnd '.cell.' lastField '{location};']);
updateGUIByGlobal(fieldname);