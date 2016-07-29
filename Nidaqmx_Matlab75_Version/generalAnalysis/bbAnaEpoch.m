function bbAnaEpoch(roi, epoch)
	
	global state
	if nargin==1
		epoch=state.epoch;
	end

	apName=['e' num2str(epoch) 'ap'];
	epspName=['e' num2str(epoch) 'epsp'];
	supraName=['e' num2str(epoch) 'supra'];
	supraCalcName=['e' num2str(epoch) 'supraCalc'];
	subName=['e' num2str(epoch) 'sub'];
	subCalcName=['e' num2str(epoch) 'subCalc'];
	
	duplicateo(['e' num2str(epoch) 'p4c1r' num2str(roi) '_avg'], apName);
	duplicateo(['e' num2str(epoch) 'p1c1r' num2str(roi) '_avg'], epspName);
	duplicateo(['e' num2str(epoch) 'p5c1r' num2str(roi) '_avg'], supraName);
	duplicateo(['e' num2str(epoch) 'p8c1r' num2str(roi) '_avg'], subName);
	
	xscale=getWave(apName, 'xscale');
	
	baselineSubtract(apName, 10, 100);
	baselineSubtract(epspName, 10, 100);
	baselineSubtract(supraName, 10, 100);
	baselineSubtract(subName, 10, 100);

	timeShiftSupra=20;		% in msec, how long after the EPSP does the AP come
	timeShiftSub=20;		% in msec, how long before the EPSP does the AP come
	
	apPeakStart=132;
	apPeakEnd=142;
	epspPeakStart=200;
	epspPeakEnd=220;
	
	apSupraPeakStart=x2pnt(apName, apPeakStart);
	apSupraPeakEnd=x2pnt(apName, apPeakEnd);
	epspSupraPeakStart=x2pnt(epspName, apPeakStart+timeShiftSupra);
	epspSupraPeakEnd=x2pnt(epspName, apPeakEnd+timeShiftSupra);
	apSubPeakStart=x2pnt(apName, epspPeakStart+timeShiftSub);
	apSubPeakEnd=x2pnt(apName, epspPeakEnd+timeShiftSub);
	epspSubPeakStart=x2pnt(epspName, epspPeakStart);
	epspSubPeakEnd=x2pnt(epspName, epspPeakEnd);
	
	apData=getWave(apName, 'data');												% gets data points from apName
	apSupraPeak=mean(apData(apSupraPeakStart:apSupraPeakEnd)); 					% defines the region of interest in apName and calculated the mean. 
	supraCalcData=getWave(epspName, 'data');									%gets data points from epspName
	epspSupraPeak=mean(supraCalcData(epspSupraPeakStart:epspSupraPeakEnd));		% defines the region of intrest in epspName and calculated the mean
	startPt=x2pnt(apName, timeShiftSupra);										%defines the start	
	endPt=length(apData)-startPt;												%defines the end
	supraCalcData(startPt:end)=supraCalcData(startPt:end)+apData(1:endPt+1); 	% adds data between the appropriate time points
	supraCalcPeak=mean(supraCalcData(epspSupraPeakStart:epspSupraPeakEnd));		% calculates the peak for superCalcData
	
	waveo(supraCalcName, supraCalcData, 'xscale', xscale);
	plot(supraName);
	append(supraCalcName);
	setPlotProps(supraName, 'color', 'red');
	setPlotProps(supraCalcName, 'color', 'red', 'linewidth', 2);
	
	apData=getWave(apName, 'data');												%this is redundant from the above section. do I need to include it?
	apSubPeak=mean(apData(apSubPeakStart:apSubPeakEnd));						%defines region of interest and calculates the mean.
	subCalcData=getWave(epspName, 'data');										%this is also redundant from the above section. do I need to include it?
	epspSubPeak=mean(subCalcData(epspSubPeakStart:epspSubPeakEnd));				%defines region of interest and calculates the mean.
	startPt=x2pnt(apName, timeShiftSub);										%defines the start
	endPt=length(apData)-startPt;												%defines the end
	subCalcData(1:endPt+1)=subCalcData(1:endPt+1)+apData(startPt:end);			% adds data between the appropriate time points
	subCalcPeak=mean(subCalcData(epspSubPeakStart:epspSubPeakEnd));				%calculate the peak for subCalcData
	
	waveo(subCalcName, subCalcData, 'xscale', xscale);
	plot(subName);
	append(subCalcName);
	setPlotProps(subName, 'color', 'green');
	setPlotProps(subCalcName, 'color', 'green', 'linewidth', 2);
	
	supraData=getWave(supraName, 'data');
	supraPeak=mean(supraData(epspSupraPeakStart:epspSupraPeakEnd));
	subData=getWave(subName, 'data');
	subPeak=mean(subData(epspSubPeakStart:epspSubPeakEnd));

	results=[epoch epspSupraPeak apSupraPeak supraPeak supraCalcPeak epspSubPeak apSubPeak subPeak subCalcPeak];
	disp('epoch	epsp        ap        supra       supraCalc   epsp       ap           sub         subCalc');
	disp(num2str(results));
	
 	if ~iswave('pairingResults')
 		waveo('pairingResults', results)
	else
		global pairingResults
		row=find(pairingResults.data(:,1)==epoch);
		if isempty(row)
			pairingResults.data(end+1, :)=results;
		else
			pairingResults.data(row(1), :)=results;
		end
	end
	if state.files.autoSave
		names={apName epspName supraName subName subCalcName subCalcName 'pairingResults'};
		for name=names
			eval(['global ' name{1}]);
			save(fullfile(state.files.savePath, name{1}), name{1});
		end
	end

