function processMCData(ai, SamplesAcquired)
	global state gh mcData

	
% 	setPhysStatusString('Processing MC Data...');
	% extract the data
    % strucure mcData is serving same purpose as acquisition waves in
    % scanimage
	[mcData.data mcData.time]=getdata(state.phys.mcAcq.mcInputDevice);
    mcData.header = state.headerString;
    mcData.data = mcData.data(:, [state.phys.mcAcq.mcChannelOrder (length(state.phys.mcAcq.mcChannelOrder) + 1):size(mcData.data, 2)]); % reorder channels...
    mcData.cycle = state.cycle;
    disp(['*** Last Cycle Position: ' num2str(state.cycle.currentCyclePosition) ' ***']);
    

    
    



% 	% do the proper scaling
% 	channelList=get(state.phys.daq.inputDevice, 'Channel');
% 	gain=ones(1,8);
	

	

% 		% store the data in a WAVE for the current acq
% 		waveo(['dataWave' num2str(channel)], eval(['state.phys.data' num2str(channel)]), ...
% 				'xscale', [0 1/state.phys.settings.inputRate]);
%         setfield(['dataWave' num2str(channel)], 'headerString', state.headerString); 
  
		
		% store the data in a wave that contains the acq #
% 		name=physTraceName(channel, state.files.lastAcquisition);
        name = [state.files.baseName 'MC_' num2str(state.files.lastAcquisition)];

% 		setfield(name, 'headerString', state.headerString);

		% auto save to disk?
		

		if state.files.autoSave
			%eval(['global ' name]);
    %        setWaveUserDataField(name, fieldName, fieldValue);0 %    This
    %        is where to add data fields for inclusion in the phys
    %        acquisition User Data....
			save(fullfile(state.files.savePath, name), 'mcData');
%			disp(['	*** Saved to disk : ' fullfile(state.files.savePath, name)]);
        end
        
        % filter data
        state.phys.mcAcq.filteredData = mcFilter(mcData.data(:, 1:state.phys.mcAcq.mcNChannels),...
            state.phys.mcAcq.globalLowPass, state.phys.mcAcq.globalHighPass, state.phys.mcAcq.mcInputRate);

        % determine Std Dev and mean of data
        state.phys.mcAcq.MUA.SD = std(state.phys.mcAcq.filteredData(:, 1:state.phys.mcAcq.mcNChannels));
        state.phys.mcAcq.MUA.means = mean(state.phys.mcAcq.filteredData(:, 1:state.phys.mcAcq.mcNChannels));
        

        
        % display data
        state.phys.mcAcq.displayXData = mcData.time;
        state.phys.mcAcq.displayData = mcData.data;
        mcAcqUpdateChannelPlots;
        % MUA analysis
        mcAcqMUA;
% 
% 		% online averaging?
% 		if getfield(state.phys.settings, ['avg' num2str(channel)])
% 			if state.cycle.useCyclePos
% 				avgName=physAvgName(state.epoch, channel, state.cycle.lastPositionUsed);
% 			else
% 				avgName=physAvgName(state.epoch, channel, state.cycle.lastPulseUsed0);
% 			end
% 			avgin(name, avgName);
%             baseline(avgName, 200); %FS MOD
%             disp('*** Averages are baselined at 200msec- this is hard coded in processPhysData ***');
% 			if state.files.autoSave
% 				eval(['global ' avgName]);
% 				save(fullfile(state.files.savePath, avgName), avgName);
% %				disp(['	*** Saved to disk : ' fullfile(state.files.savePath, avgName)]);
% 			end
% 		end% 
% 		
% 	end
	

	
	

