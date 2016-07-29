function restoreUncageInfo(posNum)
	
	global state
	
	try
		if size(state.internal.uncageInfo,2)~=19
			beep
			disp(['ERROR : restoreUncageInfo position' num2str(posNum) 'does not exist' ]);
			return
		end
	catch
		beep
		disp(['ERROR : restoreUncageInfo position' num2str(posNum) 'does not exist' ]);
		return
	end
	
	if size(state.internal.uncageInfo,1)<posNum
		beep
		disp(['ERROR : restoreUncageInfo position' num2str(posNum) 'does not exist' ]);
		return
	end
			
	state.acq.zoomFactor = state.internal.uncageInfo(posNum, 1);
	state.acq.scanRotation  = state.internal.uncageInfo(posNum, 2);
	state.acq.postRotOffsetX = state.internal.uncageInfo(posNum, 3);
	state.acq.postRotOffsetY  = state.internal.uncageInfo(posNum, 4);
	state.blaster.XList(1) = state.internal.uncageInfo(posNum, 5); 
	state.blaster.YList(1)  = state.internal.uncageInfo(posNum, 6);
	state.blaster.indexXList(1) = state.internal.uncageInfo(posNum, 7);
	state.blaster.indexYList(1) = state.internal.uncageInfo(posNum, 8);
	state.blaster.indexList(1) = state.internal.uncageInfo(posNum, 9);
	state.blaster.allConfigs{2, 2}(1, 3) = state.internal.uncageInfo(posNum, 10);
	state.blaster.allConfigs{2, 2}(1, 4) = state.internal.uncageInfo(posNum, 11);
	state.blaster.allConfigs{2, 2}(1, 5) = state.internal.uncageInfo(posNum, 12);
	state.internal.trackerX0 = state.internal.uncageInfo(posNum, 13);
	state.internal.trackerY0 = state.internal.uncageInfo(posNum, 14);
	state.acq.scanShiftX = state.internal.uncageInfo(posNum, 15);
	state.acq.scanShiftY = state.internal.uncageInfo(posNum, 16);
	state.acq.pixelShiftX = state.internal.uncageInfo(posNum, 17);
	state.acq.pixelShiftY = state.internal.uncageInfo(posNum, 18);
	state.internal.refShiftX = state.internal.uncageInfo(posNum, 19);
	state.internal.refShiftY = state.internal.uncageInfo(posNum, 20);

	state.acq.trackerReference=state.internal.trackerReferences{posNum};
	state.acq.trackerReferenceAll	= state.acq.trackerReferencesAll{posNum};

	try
		turnOffMotorButtons;
		gotoPosition(posNum);
	catch
	end
	turnOnMotorButtons;	

	updatePositionDisplay;
	updateBlasterConfigLineDisplay;

	state.internal.needNewRotatedMirrorOutput=1;
	state.internal.needNewPcellRepeatedOutput=1;
	applyChangesToOutput;
	updateReferenceImage;
	updateGUIByGlobal('state.acq.scanRotation');
	updateGUIByGlobal('state.acq.zoomFactor');
	
	disp(['*** LOADED POSITION ' num2str(posNum) ]);
	

		
function updatePositionDisplay
	global state
	state.blaster.X = state.blaster.XList(state.blaster.displayPos);
	updateGuiByGlobal('state.blaster.X');
	state.blaster.Y = state.blaster.YList(state.blaster.displayPos);
	updateGuiByGlobal('state.blaster.Y');
	
function updateBlasterConfigLineDisplay
	global state
	state.blaster.linePos=state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 1);
	state.blaster.linePat=state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 2);
	state.blaster.lineWidth=state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 3);
	state.blaster.linePower1=state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 4);
	state.blaster.linePower2=state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 5);
	state.blaster.lineTilerActive=state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 6);
	updateGuiByGlobal('state.blaster.linePos');
	updateGuiByGlobal('state.blaster.linePat');
	updateGuiByGlobal('state.blaster.lineWidth');
	updateGuiByGlobal('state.blaster.linePower1');
	updateGuiByGlobal('state.blaster.linePower2');
	updateGuiByGlobal('state.blaster.lineTilerActive');