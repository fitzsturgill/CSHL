function saveCubeAs
global state
[fname, pname] = uiputfile('*.tif', 'Choose file to save...');
if isnumeric(fname)
    return
else
    printHandleToFile(state.imageProc.cubeHandle, [pname fname]);
end
