duration = length(TE.Photometry.xData) / TE.Photometry.sampleRate;
TE = addPupilometryToTE(TE, 'duration', duration, 'zeroField', 'Cue',  'frameRate', 60, 'frameRateNew', 20);