function avgData = phAverageFromTE(TE, trials, ch, varargin)
% output - avgData with fields phMean, phSTD, phSEM, and N

% !!!!!! trials- may be cell array of trials for different conditions !!!!!!
    defaults = {...
        'PhotometryField', 'Photometry';...
        'FluorDataField', 'dFF';...
        'window', [];... % averaging window around alignment point
        'zeroTimes', [];... % empty- inferred from Photometry xData OR cell array or vector containing zero times relative to Bpod trial start
        'referenceFromEnd', 0;... % relevent when zeroTimes are supplied as a cell array (e.g. if you want to align to the beginning or end of a bpod state)        
        };    
    [s, ~] = parse_args(defaults, varargin{:});



% 
%     xData = TE.(Photometry).xData;
%     if ~isempty(s.window) && ~isempty(s.zeroTimes)
%         allTrialsZero = zeroTimes(1) - TE.(Photometry).startTime(1);
%         startP = max(1, bpX2pnt(s.window(1) + allTrialsZero, TE.(Photometry).sampleRate));
%         endP = min(length(xData), bpX2pnt(s.window(2) + allTrialsZero, TE.(Photometry).sampleRate));
%     elseif ~isempty(s.window) % use Photometry xData to infer zero point (really using first x value in Photometry xData)
%         startP = max(1, bpX2pnt(s.window(1), TE.(Photometry).sampleRate, xData(1)));
%         endP = min(length(xData), bpX2pnt(s.window(2), TE.(Photometry).sampleRate, xData(1)));        
%     else
%         startP = 1;
%         endP = length(xData);
%         s.window = [xData(1) xData(end)];
%     end

%% rewrite xData
    xData = linspace(s.window(1), s.window(2), endP - startP + 1);  % rewrite xData
    
    Avg = NaN(length(trials), length(xData));
    STD = NaN(length(trials), length(xData));
    SEM = NaN(length(trials), length(xData));
    N = NaN(length(trials), length(xData));    
    XData = NaN(length(trials), length(xData));        
    for counter = 1:length(trials)        
        currentTrials = trials{counter};
        currentData = TE.(Photometry).data(ch).(s.FluorDataField)(currentTrials, startP:endP);
        Avg(counter,:) = nanmean(currentData);
        STD(counter,:) = std(currentData, 'omitnan');
        SEM(counter,:) = std(currentData, 'omitnan') ./ sqrt(sum(~isnan(currentData), 1));
        N(counter, :) = sum(~isnan(currentData), 1);
        XData(counter, :) = xData;
    end
    
    avgData = struct(...
        'Avg', Avg,...
        'STD', STD,...
        'SEM', SEM,...
        'N', N,...
        'xData', XData...
        );

    