function figureButtonOverCallback
global state gh
% This function will grab the current position and intensity from the figure
% and display it in the imageGUI Window.
		
	% Figure out which window you ar ein
	name = get(gcf,'Name');
	lengthOfName = length(name);
	channel = str2num(name(lengthOfName));
	
	currentPoint = recordCurrentPoint(gca);
	state.internal.currentPointX = round(currentPoint(1,1));
	updateGUIByGlobal('state.internal.currentPointX');
	state.internal.currentPointY = round(currentPoint(1,2));
	updateGUIByGlobal('state.internal.currentPointY');

	try
		if state.internal.status==2	... % we are currently focusing
				|| (state.internal.status==0 && state.internal.lastTaskDone==2) ... % we just focused
				|| (state.internal.status==4 && state.internal.lastTaskDone==2) % we're looping and just focused
			if state.internal.currentPointY>0 && state.internal.currentPointX>0
				state.internal.intensity = state.acq.focusData{channel}(state.internal.currentPointY, state.internal.currentPointX);
			else
				state.internal.intensity=0;
			end
			updateGuiByGlobal('state.internal.intensity');
		else
			if state.acq.averaging == 1 
				if strcmp('Max', name(1:3))	 % Looking at a max projection
					pos=state.internal.zSliceCounter+1;
					if pos>size(state.acq.maxData{channel},3)
						pos=size(state.acq.maxData{channel},3);
					end
					state.internal.intensity = state.acq.maxData{channel}(state.internal.currentPointY, state.internal.currentPointX, pos);
				else
					pos=state.internal.zSliceCounter+1;
					if pos>size(state.acq.acquiredData{channel},3)
						pos=size(state.acq.acquiredData{channel},3);
					end
					state.internal.intensity = state.acq.acquiredData{channel}(state.internal.currentPointY, state.internal.currentPointX, pos);
				end			
			else
				if strcmp('Max', name(1:3)) % Looking at a max projection
					if state.acq.numberOfFrames==1
						pos=state.internal.zSliceCounter+1;
						if pos>size(state.acq.maxData{channel},3)
							pos=size(state.acq.maxData{channel},3);
						end
						state.internal.intensity = state.acq.maxData{channel}(state.internal.currentPointY, state.internal.currentPointX, pos);
					else
						beep;
					end
				else
					pos=max(1,state.internal.frameCounter+state.acq.numberOfFrames*state.internal.zSliceCounter);
					if pos>size(state.acq.acquiredData{channel},3)
						pos=size(state.acq.acquiredData{channel},3);
					end
					state.internal.intensity = state.acq.acquiredData{channel}(state.internal.currentPointY, state.internal.currentPointX, pos);
				end
			end
			updateGuiByGlobal('state.internal.intensity');
		end
	catch
		disp([' ERROR : figureButtonOverCallback : ' lasterr]);
	end
		
			
				
				
				
				
			