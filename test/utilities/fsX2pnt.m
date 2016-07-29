function p=mcX2pnt(x)
    global state
    
    deltaX=state.mcViewer.deltaX;
    startX=state.mcViewer.startX;
    endX=state.mcViewer.endX;
    
    p=min(max(round(1+(x-startX)/deltaX), 1), length(startX:deltaX:endX));