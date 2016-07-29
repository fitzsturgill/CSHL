function makeNewPcellPowerOutput

	global state
	state.acq.pcellPowerOutput=[];
	for pcellCounter=1:state.pcell.numberOfPcells
		powerLevel=getfield(state.pcell, ['pcellOffset' num2str(pcellCounter)]);
		offsetPower = powerToPcellVoltage(powerLevel, pcellCounter);
		if state.blaster.active & state.blaster.blankImaging
			offsetPower=powerToPcellVoltage(0, pcellCounter);
%			disp('makeNewPcellPowerOutput : blanking imaging output');
		end
		flybackPower = powerToPcellVoltage(getfield(state.pcell, ['pcellFlyBack' num2str(pcellCounter)]), pcellCounter);
		if state.acq.dualLaserMode==1
			state.acq.pcellPowerOutput(:, pcellCounter) = state.acq.pcellBinaryOutput * (offsetPower-flybackPower) + flybackPower;
		elseif state.acq.dualLaserMode==2
			if mod(pcellCounter, 2)==1
				state.acq.pcellPowerOutput(:, pcellCounter) = state.acq.pcellBinaryOutput * (offsetPower-flybackPower) + flybackPower;
			else
				state.acq.pcellPowerOutput(:, pcellCounter) = state.acq.pcellBinaryOutputComp * (offsetPower-flybackPower) + flybackPower;
			end
		end

		if powerLevel>0
			state.acq.pcellPowerOutput(:, pcellCounter + state.pcell.numberOfPcells) = 5 * state.shutter.open;
		else
			state.acq.pcellPowerOutput(:, pcellCounter + state.pcell.numberOfPcells) = 5 * state.shutter.closed;
		end				
	end
