function data = addCSVToTE(TE, varargin)

% Adds CSV data (e.g. from Bonsai) to TE.  CSV data for every session loaded into TE
% assumed to to be saved in a folder (e.g. \EyeAvg_160814\). csv data .csv
% files (e.g. EyeAvg_0.csv) are zero indexed
 
% TE: Trial Event structure containing fields trialNumber and filename
% rootPath: path containing session data and CSV session folders

    %% optional parameters, first set defaults
    defaults = {...
        'sampleRate', 60;... % integer
        'sampleRateNew', 20;... 
        'duration', 11;...  % assumes constant trial duration across TE
        'zeroField', 'Outcome';... 
        'startField', 'PreCsRecording';... % extract time stamp w.r.t. start of Bpod trial
        'rootPath', [];...
        'normMode', 'bySession';...  % [bySession, byTrial]
        'folderPrefix', 'EyeAvg_';...
        'folderSuffix', '';...
        'filePrefix', 'EyeAvg_';...
        'dataFieldName', '';... % if empty, inferred from filePrefix
        };
    
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings    
    if isempty(s.sampleRateNew)
        s.sampleRateNew = s.sampleRate;
    end        
    
    if isempty(s.rootPath)
        [~, rootPath] = uiputfile('path', 'Choose data folder containing CSV subdirectories...');
        if rootPath == 0
            return
        else
            s.rootPath = rootPath;
        end
    else
        rootPath = s.rootPath;
    end
    cd(rootPath);    
    
    if isempty(s.dataFieldName)
        if s.filePrefix(end) == '_'
            s.dataFieldName = s.filePrefix(1:end-1);
        else
            s.dataFieldName = s.filePrefix;
        end
    end
        
    
    
    %% initialize CSV structure
    dX = 1/s.sampleRateNew;
    startX = []; % to be determined
    nSamples = s.sampleRate * s.duration; 
    normFields = {...
        s.dataFieldName, [s.dataFieldName 'Norm'];...
        };
    data = struct();
    data.loadSettings = s; % scalar


    data.xData = linspace(0, s.duration - dX, s.duration * s.sampleRateNew);
    nsamplesNew = length(data.xData);
    data.startTime = NaN(length(TE.filename), 1);
    data.sampleRate = NaN(length(TE.filename), 1);
    
    for counter = 1:numel(normFields)
        data.(normFields{counter}) = NaN(length(TE.filename), nsamplesNew); % fill with NaNs
    end
    
    
    %% extract session dates in TE and generate matching CSV folder names
    sessionnames = unique(TE.filename);
    if ischar(s.folderSuffix)
        s.folderSuffix = repmat({s.folderSuffix}, length(sessionnames));
    end
    h = waitbar(0, 'Processing CSV data'); 
    
    %% assuming consistent trial durations and alignment
    startX = TE.(s.startField){1}(1) - TE.(s.zeroField){1}(1);
    data.xData = data.xData + startX;
    blStartP = 1; % just start at beginning of video for baseline
    blEndP = bpX2pnt(0, s.sampleRateNew, startX); % go to zero point        
    [p, q] = rat(s.sampleRateNew/s.sampleRate);
    for counter = 1:length(sessionnames)
        sessionname = sessionnames{counter};
        % additionally try Combined_ as a prefix 
        dataFolder = {parseFileName(sessionname, s.folderPrefix), parseFileName(sessionname, 'Combined_')}; % see subfunction
        folderFound = 0;
        for j = 1:length(dataFolder)
            dataPath = fullfile(rootPath, [dataFolder{j} s.folderSuffix{counter}], filesep); % filesep returns system file separator character
            if isdir(dataPath)
                folderFound = 1; % dataPath is correct....
                break
            end
        end
        if ~folderFound
            warning('*** CSV data directory %s does not exist ***', dataPath);
            continue
        end
        cd(dataPath);
        fs = dir([s.filePrefix '*.csv']);
        if length(fs) == 0
            warning(sprintf('*** No %s files found in %s or %s ***', s.filePrefix ,dataFolder{1} ,dataFolder{2})); 
            continue
        end
        fileList = {};
        [fileList{1:length(fs)}] = fs(:).name; % see deal documentation I think...
        [fileList,ix] = sort_nat(fileList); % alphanumeric sorting
%% verify correct numbering of pupil files        
        dmDelta = seconds(diff(datetime({fs(ix).date})));
        dmDelta = dmDelta(:);

        % find matching session indices within TE either by matching numeric suffix or total number of files 
        si = find(filterTE(TE, 'filename', sessionname));
        teDelta = diff(TE.TrialStartTimestamp(si));
        teDelta = teDelta(:);
        teDelta = circshift(teDelta, -1); % shift backward to account for Bonsai's tendency to skip the first time difference (perhaps due to Bonsai 'closing' the first movie file upon the 2nd trigger occurance rather than the end of the 11 sec duration triggeredWindow)


        numFileDifference = length(dmDelta) - length(teDelta);
        
        maxDifference = 2;

        
        % ith element of correctedIx_dm contains trial index matching the ith
        % pupil.mat file
        outlierITI = 24;        
        correctedIx_dm = 1:length(dmDelta);
        if numFileDifference < 0 % there are missing pupil files            
            startingIndex = 1; 
            for dc = 1:abs(numFileDifference)
                Rsqs = zeros(length(dmDelta) - startingIndex + 1, 1);
                for shc = startingIndex:length(dmDelta)
                    testCorrection = correctedIx_dm;
                    testCorrection(shc:end) = testCorrection(shc:end) + 1;
                    % remove outliers (hard coded here as ITIs longer than
                    % outlierITI
                    clean_dmDelta = dmDelta;
                    clean_teDelta = teDelta(testCorrection);
%                     cleanIx = (clean_dmDelta < outlierITI) & (clean_teDelta < outlierITI);          
                    cleanIx = (clean_dmDelta < clean_teDelta + outlierITI) & (clean_teDelta < 300) & (clean_teDelta > 0);
                    Rsqs(shc) = corr(clean_dmDelta(cleanIx), clean_teDelta(cleanIx));
                end
                [maxRsq, mix] = max(Rsqs);
                correctedIx_dm(mix:end) = correctedIx_dm(mix:end) + 1; % 'save' the first correction/gap prior to searching for the next
                startingIndex = mix;
            end
            dmDelta_plot = dmDelta;
            teDelta_plot = teDelta(correctedIx_dm);    
        elseif numFileDifference > 0 % there are extra pupil files
            % introduce NaNs into pupil file indices to skip over them
            startingIndex = 1; 
            for dc = 1:abs(numFileDifference)
                Rsqs = zeros(length(dmDelta) - startingIndex + 1, 1);
                for shc = startingIndex:length(dmDelta)
                    testCorrection = correctedIx_dm;
                    testCorrection(shc) = NaN;
                    clean_dmDelta = dmDelta(~isnan(testCorrection));
                    clean_dmDelta = clean_dmDelta(1:length(teDelta));
                    clean_teDelta = teDelta;
%                     cleanIx = (clean_dmDelta < outlierITI) & (clean_teDelta < outlierITI) & (clean_teDelta > 0);
                    cleanIx = (clean_dmDelta < clean_teDelta + outlierITI) & (clean_teDelta < 300) & (clean_teDelta > 0);
                    Rsqs(shc) = corr(clean_dmDelta(cleanIx), clean_teDelta(cleanIx));                   
                end
                [maxRsq, mix] = max(Rsqs);
                correctedIx_dm(mix) = NaN; % 'save' the first extra pupil file prior to searching for the next
                correctedIx_dm(min(mix+1, length(correctedIx_dm)):end) = correctedIx_dm(min(mix+1, length(correctedIx_dm)):end) - 1;
                startingIndex = mix;
            end
            dmDelta_plot = dmDelta(~isnan(correctedIx_dm));
            teDelta_plot = teDelta;            
        else
            % remove outliers (hard coded here as ITIs longer than
            % outlierITI                          
            cleanIx = (dmDelta < outlierITI) & (teDelta < outlierITI);   
            maxRsq = corr(dmDelta(cleanIx), teDelta(cleanIx));      
            dmDelta_plot = dmDelta;
            teDelta_plot = teDelta;            
        end
        % plot the corrected trial and file delta time vectors together for
        % visual inspection
        ensureFigure(sessionname, 1);
        mcLandscapeFigSetup(gcf);
        plot(dmDelta_plot, 'b'); hold on; plot(teDelta_plot, 'r'); textBox(sprintf('session %d, diff= %d, maxRsq= %.2f', counter, numel(dmDelta) - numel(teDelta), maxRsq));
        legend({'files', 'trials'});
            
        % convert correctedIx_dm from time differences back to actual
        % indices (reverse the diff)
        correctedIx_dm = [1 correctedIx_dm + 1];

        if abs(numFileDifference) > maxDifference            
            error('*** spurious or missing .csv files detected, difference of %d files for %s you may need to run  bpCleanAndVerifyBonsai ***', numFileDifference, sessionname);
        end        
        
        if any(dmDelta < 2)
            error('there are short ITI .csv files that you need to delete or deal with possibly using bpCleanAndVerifyBonsai');
        end 
        
        if maxRsq < 0.9
            warning('*** low max Rsq = %.2f for session %s  *** \n', maxRsq, sessionname);
%             continue            
        end
        
        for i = 1:length(correctedIx_dm)
            icorrected = correctedIx_dm(i);
            if isnan(icorrected)
                continue % continue if it is an extra csv file
            end                               
            tei = si(icorrected); % TE index number 
            try
                loaded = dlmread(fileList{i}, ' ');
            catch
                fprintf('*** %s didn''t exist in %s ***\n', fileList{i},  sessionname);
                continue
            end
            loaded = loaded(:,1);
            samplesAvailable = size(loaded, 1);            
            currentSamples = ceil(samplesAvailable * p/q);
            resampled = resample(loaded, p, q);
            data.sampleRate(tei, 1) = s.sampleRateNew;
            data.(s.dataFieldName)(tei,1:min(nsamplesNew, currentSamples)) = resampled(1:min(nsamplesNew, currentSamples))';   
            data.startTime(tei, 1) = TE.(s.startField){tei}(1);     
        end
            %% baseline subtract
        for i = 1:size(normFields, 1)
            rawField = normFields{i, 1};
            normField = normFields{i, 2};
            rawData = data.(rawField)(si, :);
            switch s.normMode
                case 'byTrial'
                    data.(normField)(si, :) = bsxfun(@rdivide, rawData, nanmean(rawData(:, blStartP:blEndP), 2));
                case 'bySession'
                    data.(normField)(si, :) = rawData / nanmean(nanmean(rawData(:, blStartP:blEndP), 2), 1); % divide by scalar
            end
        end
        waitbar(counter/length(sessionnames));
    end
    close(h);
end

function dataFolder = parseFileName(filename, folderPrefix)
    % extract date and generate file name matching
    % Pupil_[year][month][day] Pupil_161020 for example (october 20th)
    parts = strsplit(filename, '_');
    months = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
    m = parts{end - 2}(1:3); % Aug from Aug16
    d = parts{end - 2}(4:end); % 16 from Aug16
    y = parts{end - 1}(3:4); % 16 from 2016
    mi = strmatch(m, months);
    m2 = num2str(mi);
    if length(m2) == 1
        m2 = ['0' m2];
    end
    if length(d) == 1
        d = ['0' d];
    end
    dataFolder = [folderPrefix y m2 d];
end