%% calculate latency, jitter, and reliability of uncued Us responses

minRewardLickRate = 2;
bl_window = [-4 -3];
response_window = [0 3]; % must start at 0
maxPeakLatency = 1;  % peak must be reached in n seconds
response_window_avg = [-0.1 0.5]; % with respect to peak point in average
endTauOff = 0.8; % use up to n seconds after us delivery to fit tau off
PhotometryField = 'PhotometryHF';
na = length(DB.animals);

s2 = struct(...
    'latency', zeros(na, 2),...
    'latencyY', zeros(na, 2),... % 1/2 max value
    'avg_delta', zeros(na, 2),... % delta peak responsee from baseline
    'avg_mean', zeros(na, 2),...  % delta avg response from baseline
    'avgMax', NaN(na, 2),... % maximum (or minimum) value of response
    'tauOff', NaN(na, 2),...
    'tauObject', [],...
    'tauData', [],...
    'tauDataX', [],...
    'avgData', [],...
    'avgDataX', [],...
    'jitter_std', zeros(na, 2),...
    'jitter_avg', zeros(na, 2),...
    'jitter_tt50', [],... % cell array to hold distributions
    'jitter_delta', [],... % cell array to hold distributions
    'delta', [],... % less confusingly named cell array to hold delta distribution
    'mean', [],... % instead of a peak measurement, a mean (~area under curve)
    'Rnoise_peak', zeros(na,1),... % noise correlations for each channel pair based upon delta to peak
    'Rnoise_mean', zeros(na,1),... % noise correlations for each channel pair based upon mean measurement
    'Rnoise_bl', zeros(na, 1)... % noise correlations during baseline period for each channel pair
    );

us_pooled = struct(...
    'rew', s2,...
    'puff', s2,...
    'shock', s2,...
    'rew_cued', s2,...
    'puff_cued', s2,...
    'shock_cued', s2...
    );
us_pooled.animals = DB.animals;

trialSets = {'rew', 'puff', 'shock', 'rew_cued', 'puff_cued', 'shock_cued'};
% [rewAvgDelta, puffAvgDelta, shockAvgDelta, rewAvgDelta_cued, puffAvgDelta_cued, shockAvgDelta_cued] = deal(NaN(length(DB.animals), 2)); % hold average rew delta ZS of uncued reward trials

for counter = 1:length(DB.animals)
    animal = DB.animals{counter}
    dbLoadAnimal(DB, animal);
    Fs = TE.(PhotometryField).sampleRate;
    for tscounter = 1:length(trialSets)
        trialSet = trialSets{tscounter};
        switch trialSet
            case 'rew'
                theseTrials = find(uncuedReward & (TE.licks_us.rate > minRewardLickRate) & ismember(TE.BlockNumber, [2 3]));    
            case 'puff'
                theseTrials = find(punishTrials & uncuedTrials);   
            case 'shock'
                theseTrials = find(shockTrials & uncuedTrials);
            case 'rew_cued'
                theseTrials = find(rewardTrials & Odor2Valve1Trials & (TE.licks_us.rate > minRewardLickRate) & (TE.licks_cs.rate > minCueLickRate) & ismember(TE.BlockNumber, [2 3])); % exclude shock or extinction days
            case 'puff_cued'
                theseTrials = find(punishTrials & Odor2Valve2Trials);
            case 'shock_cued'
                theseTrials = find(shockTrials & Odor2Valve2Trials);
        end
             
    
        [tt50, deltaMax, deltaMean, bl_all] = deal(NaN(length(theseTrials), 2));
        for channel = 1:2
            [responseData, xData] = phAlignedWindow(TE, theseTrials, channel, 'zeroTimes', TE.Us, 'window', response_window, 'FluorDataField', 'ZS', 'PhotometryField', PhotometryField);
            [blData, blx] = phAlignedWindow(TE, theseTrials, channel, 'zeroTimes', TE.Us, 'window', bl_window, 'FluorDataField', 'ZS', 'PhotometryField', PhotometryField);
            bl = mean(blData, 2);
            response_avg = nanmean(responseData);
            % check if response is positive or negative, assume positive
            % for reward (so far this is correct)
            minVal = min(response_avg(1:bpX2pnt(maxPeakLatency, Fs)));
            maxVal = max(response_avg(1:bpX2pnt(maxPeakLatency, Fs)));
            if strcmp(trialSet, 'rew') || (maxVal - mean(bl)) > (mean(bl) - minVal)
                [peakVal, iAvg] = max(response_avg(1:bpX2pnt(maxPeakLatency, Fs)));
                avgDelta = peakVal - mean(bl);
                direction = 1;
                [M, I] = max(responseData, [], 2);                
                
                % fit exponential
                fitData = response_avg(iAvg:min(bpX2pnt(xData(iAvg) + endTauOff, Fs), length(response_avg)));
                fitXData = xData(iAvg:min(bpX2pnt(xData(iAvg) + endTauOff, Fs), length(response_avg)));            
                model = 'a + b * exp(-1/c*x)';
                fo = fitoptions('Method', 'NonlinearLeastSquares',...
                    'Upper', [min(response_avg), Inf, Inf],...
                    'Lower', [min(response_avg), 0, 0],...
                    'StartPoint', [min(response_avg), range(response_avg), 1]...
                    );
                ft = fittype(model, 'options', fo);
                [fitobject, gof, output] = fit(fitXData(:) - fitXData(1), fitData(:), ft, fo);
                us_pooled.(trialSet).tauOff(counter,channel) = fitobject.c;
                us_pooled.(trialSet).tauObject{counter,channel} = fitobject;
                us_pooled.(trialSet).tauData{counter,channel} = fitobject.a + fitobject.b * exp(-1/fitobject.c * (fitXData - fitXData(1)));
                us_pooled.(trialSet).tauDataX{counter,channel} = fitXData;
            else
                [peakVal, iAvg] = min(response_avg(1:bpX2pnt(maxPeakLatency, Fs)));
                avgDelta = peakVal - mean(bl);
                direction = -1;
                if strcmp(trialSet, 'rew')
                    disp([animal 'wtf']);
                end
                [M, I] = min(responseData, [], 2);                
            end

            
            us_pooled.(trialSet).avgDataX{counter,channel} = xData;
            us_pooled.(trialSet).avgData{counter,channel} = response_avg;
                                                                           
            deltaMax(:,channel) = M - bl;
            bl_all(:,channel) = bl;
            level50 = bl + deltaMax(:,channel) * 0.5;
            
            [ind, to] = crossing(response_avg, xData, mean(bl) + avgDelta/2);
            us_pooled.(trialSet).latency(counter, channel) = to(1);
            us_pooled.(trialSet).latencyY(counter, channel) = avgDelta/2;
            us_pooled.(trialSet).avg_delta(counter,channel) = avgDelta; % can be negative or positive
            
            % find time to 50% per trial and magnitude of trial response
            thresh = responseData >= level50;
            crossings = [false(size(thresh,1), 1) thresh(:,2:end) - thresh(:,1:end-1)];
            crossings = crossings > 0; % only want upward crossings
            [row, col] = find(crossings);
            pointMatrix = NaN(size(crossings));
            pointMatrix(sub2ind(size(crossings), row,col)) = col;
            valid = isfinite(min(pointMatrix,[],2)); % only include trials where crossing is found
            theseCrossings = min(pointMatrix(valid,:),[],2);
            responseDataValid = responseData(valid,:);
            x2 = xData(theseCrossings)';
            x1 = xData(theseCrossings - 1)';
            y2 = responseDataValid(sub2ind(size(responseDataValid), (1:length(theseCrossings))', theseCrossings));
            y1 = responseDataValid(sub2ind(size(responseDataValid), (1:length(theseCrossings))', theseCrossings - 1)); 
            m = (y2 - y1) ./ (x2 - x1);
            x = (level50(valid) - y1 + m .* x1) ./ m;
            tt50(valid,channel) = x;      
            
            
            [meanData, ~] = phAlignedWindow(TE, theseTrials, channel, 'zeroTimes', cellfun(@(x) x(1) + xData(iAvg), TE.Us), 'window', response_window_avg, 'FluorDataField', 'ZS', 'PhotometryField', PhotometryField);            
            peakMean = mean(meanData, 2);
            avgMean = nanmean(peakMean);
            us_pooled.(trialSet).avg_mean(counter, channel) = nanmean(peakMean - bl);
            deltaMean(:,channel) = peakMean - bl;
            
    %         ensureFigure('test', 1);
    %         trials = randperm(length(theseTrials), 16);
    %         for counter = 1:16
    %             trial = trials(counter);
    %             subplot(4,4,counter); hold on;
    %             plot(xData, responseData(trial,:), '-*');
    %             plot(blx, blData(trial,:));
    %             scatter(tt50(trial, 2), level50(trial));
    %             plot(get(gca, 'XLim'), [level50(trial) level50(trial)], '--');
    %         end
        end       
        us_pooled.(trialSet).jitter_tt50{counter} = tt50;
        us_pooled.(trialSet).jitter_std(counter, :) = nanstd(tt50);
        us_pooled.(trialSet).jitter_avg(counter, :) = nanmean(tt50);
        us_pooled.(trialSet).jitter_delta{counter} = deltaMax;
        us_pooled.(trialSet).delta{counter} = deltaMax;
        us_pooled.(trialSet).mean{counter} = deltaMean;
        us_pooled.(trialSet).Rnoise_peak(counter) = corr(deltaMax(:,1), deltaMax(:,2));
        us_pooled.(trialSet).Rnoise_mean(counter) = corr(deltaMean(:,1), deltaMean(:,2));
        us_pooled.(trialSet).Rnoise_bl(counter) = corr(bl_all(:,1), bl_all(:,2));
    end
end

save(fullfile(savepath, 'us_pooled.mat'), 'us_pooled');
disp(['*** saving: ' fullfile(savepath, 'us_pooled.mat') ' ***']);