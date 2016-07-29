function processPhysData(ai, SamplesAcquired)
	global state gh
	
	setPhysStatusString('Processing Phys Data...');
	
	% extract the data
	data=getdata(state.phys.daq.inputDevice);

	% do we need to abort?
	if state.phys.internal.abort
%		abortPhysiology;
		return
	end

	% do the proper scaling
	channelList=get(state.phys.daq.inputDevice, 'Channel');
	gain=ones(1,8);
	
	if state.phys.settings.channelType0>1
		if state.phys.settings.currentClamp0
			gain(1)=state.phys.settings.mVPerVIn0/state.phys.settings.inputGain0;
			gain(3)=state.phys.settings.pAPerVIn0;
		else
			gain(1)=state.phys.settings.pAPerVIn0/state.phys.settings.inputGain0;
			gain(3)=state.phys.settings.mVPerVIn0;
		end
	end
			
	if state.phys.settings.channelType1>1
		if state.phys.settings.currentClamp1
			gain(2)=state.phys.settings.mVPerVIn1/state.phys.settings.inputGain1;
			gain(4)=state.phys.settings.pAPerVIn1;
		else
			gain(2)=state.phys.settings.pAPerVIn1/state.phys.settings.inputGain1;
			gain(4)=state.phys.settings.mVPerVIn1;
		end
	end
	
	% record the trigger time
	eventLog=get(state.phys.daq.inputDevice, 'EventLog');
	f=find(strcmp({eventLog.Type}, 'Trigger'));
	if isempty(f)
		disp('*** processPhysData: ERROR:  No trigger information returned');
		abortPhysiology;
		return
	end
	state.phys.internal.triggerClock=eventLog(f(1)).Data.AbsTime;
	updateMinInCell(state.phys.internal.triggerClock);
	
% At this point, we got the data, scaled it and recorded the trigger time.

% Below, put it in data waves, calculate Rs if appropriate, store headerString,
% save, average, and, if desired, kill
	%state.db.wave_id=[];
	state.phys.internal.newWaves={};
	for counter=1:length(channelList)
		channel=get(channelList(counter), 'HwChannel');			% what DA #
		eval(['state.phys.data' num2str(channel) ' = gain(channel+1)*data(:,counter)'';']); % store the data
		
		if channel==0 | channel==1				% do series resistance processing for 1st 2 channels
			if eval(['state.phys.settings.channelType' num2str(channel)])>1		% if channel is a clamp
				% if channel is in V-clamp and a RS check pulse is selected
				if ~eval(['state.phys.settings.currentClamp' num2str(channel)])	& state.cycle.VCRCPulse 
					[rin, rs, cm, calcRsError]=calcrs(eval(['state.phys.data' num2str(channel)]), ...		% data
						1/state.phys.settings.inputRate, ...									% dx
						state.phys.pulses.amplitudeList(state.cycle.VCRCPulse), ...		% amp
						state.phys.pulses.delayList(state.cycle.VCRCPulse), ...			% pulse start
						state.phys.pulses.pulseWidthList(state.cycle.VCRCPulse), ...		% pulse width
						0, 0.95*state.phys.pulses.delayList(state.cycle.VCRCPulse));		% baseline start and end
					if calcRsError
						disp('processPhysData: calcRs returned an error');
					end
					eval(['state.phys.cellParams.rm' num2str(channel) '=round(10*rin)/10;']);
					eval(['state.phys.cellParams.rs' num2str(channel) '=round(10*rs)/10;']);
					eval(['state.phys.cellParams.cm' num2str(channel) '=round(10*cm)/10;']);					
					eval(['global physCellRm' num2str(counter-1)]);
					eval(['physCellRm' num2str(counter-1) '(state.files.lastAcquisition)=rin;'])
					eval(['global physCellRs' num2str(counter-1)]);
					eval(['physCellRs' num2str(counter-1) '(state.files.lastAcquisition)=rs;'])
					eval(['global physCellCm' num2str(counter-1)]);
					eval(['physCellCm' num2str(counter-1) '(state.files.lastAcquisition)=cm;'])
				else
					eval(['state.phys.cellParams.rm' num2str(channel) '=NaN;']);
					eval(['state.phys.cellParams.rs' num2str(channel) '=NaN;']);
					eval(['state.phys.cellParams.cm' num2str(channel) '=NaN;']);	
				end
				if state.files.autoSave % ADD SAVING OF PHYS PARAMETERS HERE
%                 	save(fullfile(state.files.savePath, ['physCellCm' num2str(counter-1)]), ['physCellCm' num2str(counter-1)]);
% 					save(fullfile(state.files.savePath, ['physCellRs' num2str(counter-1)]), ['physCellRs' num2str(counter-1)]);
% 					save(fullfile(state.files.savePath, ['physCellRm' num2str(counter-1)]), ['physCellRm' num2str(counter-1)]);
                end
				updateGuiByGLobal(['state.phys.cellParams.rm' num2str(channel)]);
				updateGuiByGLobal(['state.phys.cellParams.rs' num2str(channel)]);
				updateGuiByGLobal(['state.phys.cellParams.cm' num2str(channel)]);
			end
		end

		% store the data in a WAVE for the current acq
		waveo(['dataWave' num2str(channel)], eval(['state.phys.data' num2str(channel)]), ...
				'xscale', [0 1/state.phys.settings.inputRate]);
        setfield(['dataWave' num2str(channel)], 'headerString', state.headerString); 
%MOD Fitz FS- start
%         try
%             truncData=eval(['state.phys.data' num2str(channel) '(1, 1:6000)']);
%         catch
%             truncData=eval(['state.phys.data' num2str(channel)]);
%         end
%         truncData(1, end)=truncData(1, end) + 5; %scale
%         truncData(1, end-1) = truncData(1, end-1) - 20;
% 		waveo(['dataWaveTruncated' num2str(channel)], truncData, ...
% 				'xscale', [0 1/state.phys.settings.inputRate]);
% 		setfield(['dataWaveTruncated' num2str(channel)], 'headerString', state.headerString);
%MOD Fitz FS- end        
		
		% store the data in a wave that contains the acq #
		name=physTraceName(channel, state.files.lastAcquisition);
		state.phys.internal.newWaves{end+1}=name;
		waveo(name, eval(['state.phys.data' num2str(channel)]), ...
			'xscale', [0 1/state.phys.settings.inputRate]);
		setfield(name, 'headerString', state.headerString);

		% auto save to disk?
		
		eval(['global ' name]);
		
		if state.files.autoSave
			%eval(['global ' name]);
    %        setWaveUserDataField(name, fieldName, fieldValue);0 %    This
    %        is where to add data fields for inclusion in the phys
    %        acquisition User Data....
			save(fullfile(state.files.savePath, name), name);
%			disp(['	*** Saved to disk : ' fullfile(state.files.savePath, name)]);
		end

		%save to database?
% 		if state.db.conn~=0 & state.db.saveWaves 
% 			eval(['state.db.phys.wave_ids(' num2str(counter) ')=storeWaveDB(' name ');']);
% 			state.db.phys.channels(counter)=channel;
% 		end
% 		
		% online averaging?
		if getfield(state.phys.settings, ['avg' num2str(channel)])
			if state.cycle.useCyclePos
				avgName=physAvgName(state.epoch, channel, state.cycle.lastPositionUsed);
			else
				avgName=physAvgName(state.epoch, channel, state.cycle.lastPulseUsed0);
			end
			avgin(name, avgName);
%             baseline(avgName, 200); %FS MOD
%             disp('*** Averages are baselined at 200msec- this is hard coded in processPhysData ***');
			if state.files.autoSave
				eval(['global ' avgName]);
				save(fullfile(state.files.savePath, avgName), avgName);
%				disp(['	*** Saved to disk : ' fullfile(state.files.savePath, avgName)]);
			end
		end% 
        %%
        % clear waves from memory to avoid memory overrun during
        % multichannel acquisition
        % FS MOD
        if state.cycle.mcOn
            evalin('base', ['clear ' name ';']);
        end
        %
	end
	
	% Make a note in the auto notebook and save it
	addEntryToNotebook(2, ...
		[datestr(clock,13) ' (' num2str(state.phys.cellParams.minInCell0) ' min): Acq # ' num2str(state.files.lastAcquisition) ...
			' CycPos ' num2str(state.cycle.currentCyclePosition) ' Rep ' num2str(state.cycle.repeatsDone) ...
			' Patterns ' num2str(state.cycle.lastPulseUsed0) ', ' num2str(state.cycle.lastPulseUsed1) ]);
        

    
% FS MOD
    if state.cycle.mcOn
         processMCData;
    end
% END MOD

	timerRegisterPackageDone('Physiology');
	
	

