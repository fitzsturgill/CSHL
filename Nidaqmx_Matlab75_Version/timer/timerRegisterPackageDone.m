function timerRegisterPackageDone(package)
	global state gh

%	disp(['timerRegisterPackageDone: ' package]);
	timerSetPackageStatus(0, package);
	if isempty(find(state.timer.packageStatus.*state.timer.activePackages))
		state.cycle.cycleStatus=3; 	% processing

		try
			eTime=etime(state.internal.dioTriggerTime, state.phys.cellParams.breakInClock0)/60;
		catch
			eTime=etime(state.internal.dioTriggerTime, state.internal.startupTime)/60;
		end
	
		global timerAcqTime

		if ~iswave('timerAcqTime')
			wave('timerAcqTime', []);
		end
		s=size(timerAcqTime.data, 2);
		if state.files.lastAcquisition > s+1
			timerAcqTime.data(s+1:state.files.lastAcquisition)=NaN;
		end			
		timerAcqTime.data(state.files.lastAcquisition)=eTime;

		if ~isempty(state.internal.excelChannel) & state.files.autoSave
			try
				ddepoke(state.internal.excelChannel, ['r' num2str(25 + state.files.lastAcquisition) 'c1:r' num2str(25 + state.files.lastAcquisition) 'c5'], ...
					[...
						state.files.lastAcquisition ...
						state.internal.triggerTimeInSeconds...
						state.epoch...
						0 ...
						state.cycle.currentCyclePosition...
					]);
				ddepoke(state.internal.excelChannel, ['r' num2str(25 + state.files.lastAcquisition) 'c4'], state.cycle.cycleName);
				ddepoke(state.internal.excelChannel, ['r' num2str(25 + state.files.lastAcquisition) 'c6'], state.analysis.setupName);
			catch
				disp('timerRegisterPackageDone : Could not link to excel');
			end
		end
		try
            if state.db.conn~=0
                addRecordByTable('acquisition'); %set state.db.acq_id to serial of      
                state.db.acq_id=getLastSerialInsert('acquisition', 'acq_id');
            end
        catch
        end
		timerCallPackageFunctions('Process');
		
		if state.files.autoSave		% BSMOD - Check status of autoSave option
%             if state.internal.saveHeaderAsTxt   %TN
%                 saveHeadertoTxt; %TN
%             end                  %TN
			if state.notebook.autoSaveNotes
				saveNotebooks;
            end            
            state.files.fileCounter=state.files.fileCounter+1;
			updateGUIByGlobal('state.files.fileCounter');
			updateFullFileName;
		end

		state.internal.lastTimeDelay=state.cycle.nextTimeDelay;

		if state.cycle.loopingStatus==0		% not a loop, just a single 
			set(gh.timerMainControls.doOne, 'String', 'DO ONE');
			set([gh.timerMainControls.doOne gh.timerMainControls.loop], 'Visible', 'on');
			state.cycle.cycleStatus=0;
        else
            disp('resuming Loop');
			timerResumeLoop;
			state.cycle.cycleStatus=0;
		end
	end
	
	