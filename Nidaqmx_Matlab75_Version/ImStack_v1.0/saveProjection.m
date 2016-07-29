function saveProjection
global state gh
% Saves the projection or averaging to disk...

filename = get(gh.maxProjectionGUI.fileName, 'String');
value = get(gh.maxProjectionGUI.fileName, 'Value');
if ischar(filename)
    [path,name,ext,ver] = fileparts(filename);
    if isdir(path)
        cd(path);
    end
else
    filename=filename{value};
    [path,name,ext,ver] = fileparts(filename);
    if isdir(path)
        cd(path);
    end
end

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
    maxProjection = collapse(state.imageProc.cell.currentImage{value}(:,:, ...
        state.imageProc.maxStart:state.imageProc.maxEnd), direction,'average');
else
    maxProjection = collapse(state.imageProc.cell.currentImage{value}(:,:, ...
        state.imageProc.maxStart:state.imageProc.maxEnd), direction);
end
%  write the file to disk
[fname, pname] = uiputfile('*.tif', ['Save Projection of ' filename ' along the ' direction  ...
        ' as...' ]);
newfilename = [pname fname];
arrayToTiff(maxProjection, newfilename, '');
