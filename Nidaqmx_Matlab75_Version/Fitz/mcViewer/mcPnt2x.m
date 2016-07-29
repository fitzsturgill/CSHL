function x=mcPnt2x(p)
    global state
    
    deltaX=state.mcViewer.deltaX;
    startX=state.mcViewer.startX;
    
    x = startX + deltaX * (p - 1);
	