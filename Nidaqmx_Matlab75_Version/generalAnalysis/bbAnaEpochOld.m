function bbAnaEpoch(roi, epoch)
	
	if nargin==1
		global state
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

	timeShiftSupra=10;		% in msec, how long after the EPSP does the AP come
	timeShiftSub=20;		% in msec, how long before the EPSP does the AP come
	
	peakStart=x2pnt(apName, 152);
	peakEnd=x2pnt(apName, 162);
	
	apData=getWave(apName, 'data');
	apPeak=mean(apData(peakStart:peakEnd));
	supraCalcData=getWave(epspName, 'data');
	epspPeak=mean(supraCalcData(peakStart:peakEnd));
	startPt=x2pnt(apName, timeShiftSupra);
	endPt=length(apData)-startPt;
	supraCalcData(startPt:end)=supraCalcData(startPt:end)+apData(1:endPt+1);
	supraCalcPeak=mean(supraCalcData(peakStart:peakEnd));

	waveo(supraCalcName, supraCalcData, 'xscale', xscale);
	plot(supraName);
	append(supraCalcName);
	setPlotProps(supraName, 'color', 'red');
	setPlotProps(supraCalcName, 'color', 'red', 'linewidth', 2);
	
	subCalcData=getWave(epspName, 'data');
	startPt=x2pnt(apName, timeShiftSub);
	endPt=length(apData)-startPt;
	subCalcData(1:endPt+1)=subCalcData(1:endPt+1)+apData(startPt:end);
	subCalcPeak=mean(subCalcData(peakStart:peakEnd));
	
	waveo(subCalcName, subCalcData, 'xscale', xscale);
	plot(subName);
	append(subCalcName);
	setPlotProps(subName, 'color', 'green');
	setPlotProps(subCalcName, 'color', 'green', 'linewidth', 2);

	supraData=getWave(supraName, 'data');
	supraPeak=mean(supraData(peakStart:peakEnd));
	subData=getWave(subName, 'data');
	subPeak=mean(subData(peakStart:peakEnd));
	
	results=[epoch epspPeak apPeak supraPeak supraCalcPeak subPeak subCalcPeak];
	disp('epoch    epsp       ap         supra       supraCalc   sub         subCalc');
	disp(num2str(results))
	
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

