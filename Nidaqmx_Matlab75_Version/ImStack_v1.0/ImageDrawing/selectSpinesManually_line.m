function selectSpinesManually_line
	global gh state

	peaks=[];

	figure(gh.spineGUI.mainFigure);
%	data=improfile';
	data=smooth2(improfile',3);
	waveo('FWHMdata', data);

	offset=min(data);		% use this line to use a local minimum as the offset
							% use the next line to use the true PMT offset
%	offset = valueFromHeaderString('state.acq.pmtOffsetChannel1', state.imageProc.spine.headerTop) * valueFromHeaderString('state.acq.binFactor', state.imageProc.spine.headerTop);

	waveo('offsetWave', repmat(offset, 1, length(data))');

	[roix, roiy]=findPeaks(data, 2, offset, 0.3);	% 0.3 = 30% of the peak.  Change as necessary
	
	if (size(roix, 1) < 2) | (size(roix, 2) ~= 2)
		founderr=1;
		beep;
		disp('*** COULD NOT FIND 2 ROIs');
		return
	end

	for counter=1:2
		if counter<=size(roix, 1)
   			waveo(['roi_' num2str(counter) 'x'], roix(counter, :)-1);
		    waveo(['roi_' num2str(counter) 'y'], roiy(counter, :));
		else
            if iswave(['roi_' num2str(counter) 'x'])
    			setWave(['roi_' num2str(counter) 'x'], 'data', []);
            end                
            if iswave(['roi_' num2str(counter) 'y'])
    			setWave(['roi_' num2str(counter) 'y'], 'data', []);
            end                
		end			
	end
	
	hw1=abs(roix(1,2)-roix(1,1));
	[p1, i1]=max(data(roix(1,1):roix(1,2)));
	hw2=abs(roix(2,2)-roix(2,1));
	[p2, i2]=max(data(roix(2,1):roix(2,2)));
	d=abs((i1+roix(1,1)-1) - (i2+roix(2,1)-1));
	
	state.imageProc.spine.numberOfSpines = state.imageProc.spine.numberOfSpines + 5;
	updateGUIByGlobal('state.imageProc.spine.numberOfSpines');
	state.imageProc.spine.spineLengths(end+1) = hw1;
	state.imageProc.spine.spineLengths(end+1) = p1;
	state.imageProc.spine.spineLengths(end+1) = hw2;
	state.imageProc.spine.spineLengths(end+1) = p2;
	state.imageProc.spine.spineLengths(end+1) = d;

	saveSpineDataToExcel;