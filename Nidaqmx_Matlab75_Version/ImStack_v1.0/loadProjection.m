function loadProjection
global state gh

value = get(gh.maxProjectionGUI.fileName, 'Value');

if get(gh.maxProjectionGUI.XYMax, 'Value')
    direction = get(gh.maxProjectionGUI.XYMax, 'String');
elseif get(gh.maxProjectionGUI.XZMax, 'Value')
    direction = get(gh.maxProjectionGUI.XZMax, 'String');
elseif get(gh.maxProjectionGUI.YZMax, 'Value')
    direction = get(gh.maxProjectionGUI.YZMax, 'String');
else
    display('No Direction Specified: Check state.imageProc.internal.');
end
if state.imageProc.averageNotProject    %average not project
    eval(['state.imageProc.internal.maxProjection' num2str(state.imageProc.internal.maxCounter) ...
            ' = collapse(state.imageProc.cell.currentImage{value}(:,:,state.imageProc.maxStart:state.imageProc.maxEnd), direction,''average'');']);
else
    eval(['state.imageProc.internal.maxProjection' num2str(state.imageProc.internal.maxCounter) ...
            ' = collapse(state.imageProc.cell.currentImage{value}(:,:,state.imageProc.maxStart:state.imageProc.maxEnd), direction);']);
    
end
% set the first current image array to the max image created.
loadImageFromArray(['state.imageProc.internal.maxProjection' num2str(state.imageProc.internal.maxCounter)]);
state.imageProc.internal.maxCounter = state.imageProc.internal.maxCounter+1;