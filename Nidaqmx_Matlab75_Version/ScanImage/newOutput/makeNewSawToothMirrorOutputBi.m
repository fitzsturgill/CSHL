function makeNewSawToothMirrorOutputBi
	global state

	actualOutputRate = get(state.daq.grabOutput, 'SampleRate');
	state.internal.lengthOfXData = actualOutputRate*state.acq.msPerLine;
	state.internal.lineDelay = .001*state.acq.lineDelay/state.acq.msPerLine;
	state.internal.flybackDecimal= 1-state.acq.fillFraction-state.internal.lineDelay;
	
	state.internal.startOutputColumnInLine=round(state.internal.lengthOfXData*(state.internal.flybackDecimal/2))+1;;
	state.internal.endOutputColumnInLine=round(state.internal.lengthOfXData*(1-state.internal.flybackDecimal/2))+1;
	
	state.internal.startOutputFractionInLine=(state.internal.startOutputColumnInLine-1)/state.internal.lengthOfXData;
	state.internal.endOutputFractionInLine=(state.internal.endOutputColumnInLine-1)/state.internal.lengthOfXData;
	
	oneSinWave=cos(2*pi*linspace(0, 1, 2*state.internal.lengthOfXData));
	x0=round(state.internal.lengthOfXData*state.internal.flybackDecimal/2);
	x1=round(state.internal.lengthOfXData-state.internal.lengthOfXData*state.internal.flybackDecimal/2);
	oneSinWaveLin=oneSinWave;
	oneSinWaveLin(x0:x1)=linspace(oneSinWave(x0), oneSinWave(x1), x1-x0+1);
	state.internal.start1=x0;
	state.internal.end1=x1;
	
	x0=round(state.internal.lengthOfXData+state.internal.lengthOfXData*state.internal.flybackDecimal/2);
	x1=round(2*state.internal.lengthOfXData-state.internal.lengthOfXData*state.internal.flybackDecimal/2);
	oneSinWaveLin(x0:x1)=linspace(oneSinWave(x0), oneSinWave(x1), x1-x0+1);
	state.internal.start2=x0;
	state.internal.end2=x1;
	
 	state.acq.rawSawtoothMirrorOutput=repmat(oneSinWaveLin, 1, state.acq.linesPerFrame/2)';
	
	y=reshape(repmat(linspace(-1,1,state.acq.linesPerFrame),state.internal.lengthOfXData,1),1,state.internal.lengthOfXData*state.acq.linesPerFrame);
	
	y2=[y(state.internal.start1 : state.internal.lengthOfXData*state.acq.linesPerFrame+1-state.internal.start1) repmat(-1, 1, 2*state.internal.start1-2)];
	state.acq.rawSawtoothMirrorOutput(:,2)=y2';
		
	actualInputRate = get(state.daq.grabInput, 'SampleRate');
	
 	state.internal.startDataColumnInLine = ...
 		1 + round((state.internal.flybackDecimal/2+0.001*(state.acq.lineDelay+state.acq.mirrorLag)/state.acq.msPerLine)*state.internal.samplesPerLine);
	state.internal.endDataColumnInLine = state.internal.startDataColumnInLine + (state.acq.samplesAcquiredPerLine-1);	
	
	