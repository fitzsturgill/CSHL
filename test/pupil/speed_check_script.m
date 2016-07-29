tic
for counter = 1:state.pupil.nFrames
    state.pupil.currentFrame = counter;
    updateGUIByGlobal('state.pupil.currentFrame');
    pupProcessFrame;
    pupUpdateFrameFigure;
    drawnow;
end
toc
%%
tic

for counter = 1:state.pupil.nFrames
    pupFlipFrame(counter, 1);
    
end
toc
%%
tic
for counter = 1:state.pupil.nFrames
    pupFlipFrame(counter, 0);
end
toc


%%

% quick save script for lab meeting
pupilData = struct(...
    'eyeArea', NaN(1,state.pupil.nFrames),...
    'eyeMinorAxisLength', NaN(1, state.pupil.nFrames),...
    'pupDiameter', NaN(1, state.pupil.nFrames),...
    'pupCircResidual', NaN(1, state.pupil.nFrames)...
    );
ensureFigure('pupilData', 1);
dummyData = zeros(1081, 1);
ax=zeros(1,4);
ln = zeros(1,4);
ax(1) = subplot(2,2,1); ln(1) = plot(dummyData); ylabel('eye area');
ax(2) = subplot(2,2,2); ln(2) = plot(dummyData); ylabel('minor axis length');
ax(3) = subplot(2,2,3); ln(3) = plot(dummyData); ylabel('diameter');
ax(4) = subplot(2,2,4); ln(4) = plot(dummyData); ylabel('pupil radius fit residuals');
for counter = 1:state.pupil.nFrames
    pupFlipFrame(counter, 1);
    pupilData.eyeArea(counter) = state.pupil.eye.area;
    set(ln(1), 'YData', pupilData.eyeArea);
    pupilData.eyeMinorAxisLength(counter) = state.pupil.eye.minorAxisLength;
    set(ln(2), 'YData', pupilData.eyeMinorAxisLength);
    pupilData.pupDiameter(counter) = state.pupil.pupil.diameter;
    set(ln(3), 'YData', pupilData.pupDiameter);    
    pupilData.pupCircResidual(counter) = state.pupil.pupil.circResidual;    
    set(ln(4), 'YData', pupilData.pupCircResidual);    
    drawnow;
end
toc
%%
ensureFigure('pupilData', 1);
subplot(2,2,1); plot(pupilData.eyeArea);
subplot(2,2,2); plot(pupilData.eyeMinorAxisLength);
subplot(2,2,3); plot(pupilData.pupDiameter);
subplot(2,2,4); plot(pupilData.pupCircResidual);
