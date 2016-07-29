function saveOverlayAs
global gh state

try
    [path,name,ext] = fileparts(state.imageProc.fileName);
    newfilename = [path name];
    cd(path);
end

[fname, pname] = uiputfile('*.tif','Save Overlayed Image as...' );
if isnumeric(fname)
    return
else
    filename = [pname fname];
    if ishandle(state.imageProc.overlay.figure)
        printHandleToFile(state.imageProc.overlay.figure, filename);
    else
        error('saveOverlayAs: Must have an overlay displayed to save');
    end
end



