function resetAMPKmorph(num)
    global state

    state.files.fileCounter=1;
    updateGUIbyGlobal('state.files.fileCounter');
    state.acq.zStepSize=1;
    updateGUIbyGlobal('state.acq.zStepSize');
    state.acq.numberOfFrames=2;
    updateGUIbyGlobal('state.acq.numberOfFrames');
    FSSetPath_AMPK(num);
end