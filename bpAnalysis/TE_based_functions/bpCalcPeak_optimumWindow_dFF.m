function optimizePeak = bpCalcPeak_optimumWindow_dFF(TE, maxTrials, minTrials, varargin)
    % calculates peak fluorescence using a time subwindow that maximizes
    % the discriminability of fluorescence values for maxTrials (big
    % fluorescence response expected) vs. minTrials (small fluroescence
    % response expected)
    % Uses a grid search algorithm
    % syntax differs from bpCalcPeak_dFF (uses phAlignedWindow which
    % operates on TE rather than Photometry substructure)
    
    defaults = {...
        'plot', false;... % 'mean' or 'max', 'min', 'percentile', 'center' center =  center of mass normalized to window (i.e. 0.5 = middle of window)
        'FluorDataField', 'ZS';...
        'PhotometryField', 'Photometry';...
        'window', [];... % not optional
        'zeroTimes', [];...  % not optional
        'channels', [];... % if empty, use all channels present
        'minWindow', [];... % if provided, specifies a minimum width for the optimized window
        'maxWindow', [];... % if provided, specifies a maximum width for the optimized window
        'fig', [];...
        };

    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    
    assert(~isempty(s.window), 'window is required parameter');
    assert(~isempty(s.zeroTimes), 'zeroTimes is required parameter');
    assert(isfield(TE, s.PhotometryField), [s.PhotometryField ' field does not exist']);
    
    if isempty(s.minWindow)
        s.minWindow = 0;
    end
    if isempty(s.maxWindow)
        s.maxWindow = Inf;
    end
    
    if s.plot
        if isempty(s.fig)
            saveName = 'OptimizeWindow';                
            fh = ensureFigure(saveName, 1);
            mcLandscapeFigSetup(gcf);
        else
            figure(s.fig); % bring it to front
            mcLandscapeFigSetup(s.fig); % and size it
        end
    end
    
    Fs = TE.(s.PhotometryField).sampleRate;
    if isempty(s.channels)
        s.channels = TE.(s.PhotometryField).settings.channels;
    end
    commonWindow = [max(s.window(:,1)) min(s.window(:,2))]; % use the smallest common window

    
    % collect phWindow for commonCuewindow for each channel      
    nChannels = length(s.channels);
    for channel = s.channels          
    %     find peak cs response
        avgData = phAverageFromTE(TE, maxTrials, channel, 'zeroTimes', s.zeroTimes, 'window', commonWindow, 'FluorDataField', s.FluorDataField, 'PhotometryField', s.PhotometryField); 
        [mv, mix] = max(avgData.Avg);
        % if there is no peak, use point 2/3 way into window
        if (mix > 0.9 * length(avgData.Avg))
            mix = round(1/3 * length(avgData.Avg));
            mv = avgData.Avg(mix);
        end
        prePoints = (1:2:mix)';
%         postPoints = mix+1:2:min(bpX2pnt(commonWindow(2), Fs, 0), diff(commonWindow) * Fs);
        postPoints = mix+1:2:(diff(commonWindow) * Fs);
        % make matrix of pre, post points, and auROC, scaled for each pairing
        preMatrix = repmat(prePoints, 1, length(postPoints));
        postMatrix = repmat(postPoints, length(prePoints), 1);
        window_auROC = zeros(length(prePoints), length(postPoints));
        window_delta = (postMatrix - preMatrix + 1) ./ Fs;
        [phData_plus, ~] = phAlignedWindow(TE, maxTrials, channel, 'zeroTimes', s.zeroTimes, 'window', commonWindow, 'FluorDataField', s.FluorDataField, 'PhotometryField', s.PhotometryField); 
        [phData_minus, ~] = phAlignedWindow(TE, minTrials, channel, 'zeroTimes', s.zeroTimes, 'window', commonWindow, 'FluorDataField', s.FluorDataField, 'PhotometryField', s.PhotometryField); 
        for pcounter = 1:numel(preMatrix)
            [D, P] = rocarea(nanmean(phData_plus(:,preMatrix(pcounter):postMatrix(pcounter)), 2), nanmean(phData_minus(:,preMatrix(pcounter):postMatrix(pcounter)), 2), 'scale');
            window_auROC(pcounter) = D;
        end

        [m, mwix] = max(window_auROC(:) .* (window_delta(:) >= s.minWindow) .* (window_delta(:) <= s.maxWindow));

        WindowOptIx = [preMatrix(mwix) postMatrix(mwix)];
        WindowOpt = [bpPnt2x(WindowOptIx(1), Fs) bpPnt2x(WindowOptIx(2), Fs)];
        optimizePeak(channel) = bpCalcPeak_dFF(TE.(s.PhotometryField), channel, WindowOpt, TE.Cue, 'method', 'mean', 'phField', 'ZS');
        optimizePeak(channel).settings.optimumWindowIx = WindowOptIx; % in points from Cue onset
        optimizePeak(channel).settings.optimumWindow = WindowOpt; % in seconds
        optimizePeak(channel).settings.optimumWindow_avg = avgData.Avg;
        optimizePeak(channel).settings.optimumWindow_avgX = avgData.xData;
        optimizePeak(channel).settings.optimumWindow_watershedPoint = mix;  
        optimizePeak(channel).settings.optimumWindow_StartPointMatrix = avgData.xData(preMatrix);
        optimizePeak(channel).settings.optimumWindow_StopPointMatrix = avgData.xData(postMatrix);
        optimizePeak(channel).settings.optimumWindow_auROC = window_auROC;
        optimizePeak(channel).settings.optimumWindow_auROCMax = window_auROC(mwix);
        
        if s.plot
            subplot(2,nChannels,channel);  surf(optimizePeak(channel).settings.optimumWindow_StartPointMatrix,...
                optimizePeak(channel).settings.optimumWindow_StopPointMatrix,...
                window_auROC); hold on; 
            plot3(avgData.xData(preMatrix(mwix)), avgData.xData(postMatrix(mwix)), window_auROC(mwix), 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
            title(sprintf('channel %u', channel), 'Interpreter', 'none');        
            xlabel('time from cue (s)'); ylabel('time from cue (s)');
            subplot(2,nChannels,channel + nChannels);
            plot(avgData.xData, avgData.Avg); hold on;
            stem(avgData.xData(mix), avgData.Avg(mix), 'k');
            ylabel(s.FluorDataField); xlabel('time from cue (s)');
            addStimulusPatch(gca, [WindowOpt 0 0.75], sprintf('auROC=%.2f, Nmax=%u, Nmin=%u', window_auROC(mwix),...
                sum(maxTrials), sum(minTrials)), [1 0 0]);
        end        
    end