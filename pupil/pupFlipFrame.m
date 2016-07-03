function pupFlipFrame(frame, show)
    global state
    
    if nargin < 1
        frame = state.pupil.currentFrame;
    else
        state.pupil.currentFrame = frame;  
    end
    
    if nargin < 2
        show = 1; % show to actually update figures, might want to avoid updating for speed during batch processing
    end
    state.pupil.currentFrame;
    pupProcessFrame(frame);
    if show
        updateGUIByGlobal('state.pupil.currentFrame');
        pupUpdateFrameFigure;
    end