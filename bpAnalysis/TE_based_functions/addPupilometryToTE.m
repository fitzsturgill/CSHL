function TE = addPupilometryToTE(TE, varargin)

% Adds pupil data to TE.  Pupil data for every session loaded into TE
% assumed to to be saved in a folder (e.g. \Pupil_160814\). Pupil data .mat
% files (e.g. Pupil_0.mat) are zero indexed

% TE: Trial Event structure containing fields trialNumber and filename
% rootPath: path containing session data and pupil session folders

    %% optional parameters, first set defaults
    defaults = {...
        'frameRate', 60;... % integer
        'frameRateNew', 20;... 
        'duration', 11;...  % assumes constant trial duration across TE
        'zeroField', 'Us';... 
        'startField', 'PreCsRecording';... % extract time stamp w.r.t. start of Bpod trial
        'rootPath', [];...
        'normMode', 'bySession';...  % [bySession, byTrial]
        'folderSuffix', '';...
        'fillNaNs', 0;...
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings    
    if isempty(s.frameRateNew)
        s.frameRateNew = s.frameRate;
    end        
    [p, q] = rat(s.frameRateNew/s.frameRate);   % coefficients for resampling
    if isempty(s.rootPath)
        [~, rootPath] = uiputfile('path', 'Choose data folder containing pupil subdirectories...');
        if rootPath == 0
            return
        else
            s.rootPath = rootPath;
        end
    else
        rootPath = s.rootPath;
    end
    cd(rootPath);
    maxDiameter = 200; % sometimes the diameter calculation blows up if this happens, set it to NaN
    %% initialize pupil structure
    dX = 1/s.frameRateNew;
    startX = []; % to be determined
    nFrames = s.frameRate * s.duration; 
    normFields = {...
        'eyeArea', 'eyeAreaNorm';...
        'pupArea', 'pupAreaNorm';...
        'pupDiameter', 'pupDiameterNorm';...
        'pupResidual', 'pupResidualNorm';...
        };
    pupil = struct();
    pupil.loadSettings = s; % scalar
    pupil.settings = cell(length(TE.filename), 1);
%     pupil.xData = 0:dX:(nFrames - 1) * dX;
    pupil.xData = linspace(0, s.duration - dX, s.duration * s.frameRateNew);
    nFramesNew = length(pupil.xData);
    pupil.startTime = NaN(length(TE.filename), 1);
    pupil.frameRate = NaN(length(TE.filename), 1);
    pupil.blinkDetected = NaN(length(TE.filename), nFramesNew);
    
    for counter = 1:numel(normFields)
        pupil.(normFields{counter}) = NaN(length(TE.filename), nFramesNew); % fill with NaNs
    end
    
    %% extract session dates in TE and generate matching pupil folder names
    sessionnames = unique(TE.filename);
    if ischar(s.folderSuffix)
        s.folderSuffix = repmat({s.folderSuffix}, length(sessionnames));
    end
    

    %% assuming consistent trial durations and alignment
    startX = TE.(s.startField){1}(1) - TE.(s.zeroField){1}(1);
    pupil.xData = pupil.xData + startX;
    blStartP = 1; % just start at beginning of video for baseline
    blEndP = bpX2pnt(0, s.frameRateNew, startX); % go to zero point        
    
    for counter = 1:length(sessionnames)
        sessionname = sessionnames{counter};
        % try Pupil_ and Combined_ as a prefixes
        pupilFolder = {parseFileName(sessionname, 'Pupil_'), parseFileName(sessionname, 'Combined_')}; % see subfunction
        folderFound = 0;
        for j = 1:length(pupilFolder)
            pupilPath = fullfile(rootPath, [pupilFolder{j} s.folderSuffix{counter}], filesep); % filesep returns system file separator character
            if isdir(pupilPath)
                folderFound = 1; % pupilPath is correct....
                break
            end
        end
        if ~folderFound
            warning('*** Pupil directory does not exist ***');
            continue
        end
        cd(pupilPath);

        fs = dir('*upil_*.avi'); % use videos to get time stamps, then generate .mat file name later on which will be derived from the timstamped video file
        if length(fs) == 0
            warning(['*** No Pupil_ files found in ' pupilFolder]); % you need to analyze it first!
            continue
        end
        fileList = {fs(:).name};        
        [fileList,ix] = sort_nat(fileList); % alphanumeric sorting
%% verify correct numbering of pupil files        
        dmDelta = seconds(diff(datetime({fs(ix).date})));
        dmDelta = dmDelta(:);
        if any(dmDelta < 2)
            error('there are spurious pupil files that you need to delete or deal with');
        end

        % find matching session indices within TE either by matching numeric suffix or total number of files 
        si = find(filterTE(TE, 'filename', sessionname));
        teDelta = diff(TE.TrialStartTimestamp(si));
        teDelta = teDelta(:);
        teDelta = circshift(teDelta, -1); % shift backward because Bonsai apparently time-stamps the files upon the second trigger rather the end of the 1st trigger's window (e.g. 11sec post-first trigger for reversal experiment)


        numFileDifference = length(dmDelta) - length(teDelta);
        
        maxDifference = 2;
        if abs(numFileDifference) > maxDifference            
            error(['*** spurious or missing pupil.mat files detected for ' sessionname ' you may need to run  bpCleanAndVerifyBonsai ***']);
        end        
        
        % ith element of correctedIx contains trial index matching the ith
        % pupil.mat file
        correctedIx = 1:length(dmDelta);
        if numFileDifference < 0 % there are missing pupil files            
            startingIndex = 1; 
            for dc = 1:abs(numFileDifference)
                Rsqs = zeros(length(dmDelta) - startingIndex + 1, 1);
                for shc = startingIndex:length(dmDelta)
                    testCorrection = correctedIx;
                    testCorrection(shc:end) = testCorrection(shc:end) + 1;
                    Rsqs(shc) = corr(dmDelta, teDelta(testCorrection));
                end
                [maxRsq, mix] = max(Rsqs);
                correctedIx(mix:end) = correctedIx(mix:end) + 1; % 'save' the first correction/gap prior to searching for the next
                startingIndex = mix;
            end
            ensureFigure(sessionname, 1);
            mcLandscapeFigSetup(gcf);
            plot(dmDelta, 'b'); hold on; plot(teDelta(correctedIx), 'r'); textBox(sprintf('session %d, diff= %d, maxRsq= %.2f', counter, numel(dmDelta) - numel(teDelta), corr(dmDelta, teDelta(correctedIx))));
            legend({'files', 'trials'});
        else
            if numFileDifference > 0
                % RIGHT NOW I'M ASSUMING THAT EXTRA PUPIL FILES OCCUR AT THE
                % VERY END
                correctedIx = correctedIx(1:length(teDelta));
            end
            maxRsq = corr(dmDelta(1:length(teDelta)), teDelta);
            ensureFigure(sessionname, 1);
            mcLandscapeFigSetup(gcf);
            plot(dmDelta, 'b'); hold on; plot(teDelta, 'r'); textBox(sprintf('session %d, diff= %d, maxRsq= %.2f', counter, numel(dmDelta) - numel(teDelta), maxRsq));
            legend({'files', 'trials'});            
        end
        % convert correctedIx from time differences back to actual
        % indices (reverse the diff)
        correctedIx = [1 correctedIx + 1];
        
        if maxRsq < 0.9
            warning('*** low max Rsq = %.2f for session %s  *** \n', maxRsq, sessionname);
%             continue            
        end
        
        %% loop through pupil .mat files from a given session
        wb = waitbar(0, num2str(counter));
        for i = 1:length(correctedIx)
            icorrected = correctedIx(i);
            tei = si(icorrected); % TE index number               
            try
                [~,fname,~] = fileparts(fileList{i});                 
                loaded = load([fname '.mat']);
            catch
                fprintf('*** %s didn''t exist in %s ***\n', [fname '.mat'], pupilPath)
                continue;
            end
            
            %% interpolate/extrapolate NaNs if desired
            if s.fillNaNs
                fillFields = {'diameter', 'area'};
                for fcounter = 1:length(fillFields)
                    try
                        loaded.pupilData.pupil.(fillFields{fcounter}) = inpaint_nans(loaded.pupilData.pupil.(fillFields{fcounter}));
                    catch
                        warning('inpaint_nans failed in addPuupilometryToTE for %', [fname '.mat']);
                    end
                end
            end
%%          just use pre-allocated NaNs as pad
            framesToLoad = min(nFrames, loaded.pupilData.currentFrame(end)); % bonsai seems most often to add one extra frame, still handling case where a few frames are dropped, see below

            newFrames = ceil(framesToLoad * p/q); % see resample documentation, weird brackets mean 'ceil'
            pupil.frameRate(tei,1) = s.frameRateNew;
            pupil.settings{tei} = loaded.pupilData.settings;          

            pupil.eyeArea(tei, 1:newFrames) = changeRate(loaded.pupilData.eye.area(1:framesToLoad), s, p, q);
            pupil.blinkDetected(tei, 1:newFrames) = changeRate(loaded.pupilData.eye.blinkDetected(1:framesToLoad), s, p, q);
            pupil.pupArea(tei, 1:newFrames) = changeRate(loaded.pupilData.pupil.area(1:framesToLoad), s, p, q);
            pupil.pupDiameter(tei, 1:newFrames) = changeRate(loaded.pupilData.pupil.diameter(1:framesToLoad), s, p, q);
            pupil.pupResidual(tei, 1:newFrames) = changeRate(loaded.pupilData.pupil.circResidual(1:framesToLoad), s, p, q);
            pupil.startTime(tei, 1) = TE.(s.startField){tei}(1);            
            waitbar(i/length(si));
        end
        close(wb);
        %% remove outliers
        pupil.pupDiameter(pupil.pupDiameter > maxDiameter) = NaN;
        
        %% baseline subtract
        for i = 1:size(normFields, 1)
            rawField = normFields{i, 1};
            normField = normFields{i, 2};
            rawData = pupil.(rawField)(si, :);
            switch s.normMode
                case 'byTrial'
                    pupil.(normField)(si, :) = bsxfun(@rdivide, rawData, nanmean(rawData(:, blStartP:blEndP), 2));
                case 'bySession'
                    pupil.(normField)(si, :) = rawData / nanmean(nanmean(rawData(:, blStartP:blEndP), 2), 1); % divide by scalar
            end
        end
                    
            
    end
    TE.pupil = pupil;

% save(fullfile(rootPath, 'WTF.mat'), 'wtf');
end

    %% subfunctions
function out = changeRate(data, s, p, q)
    if s.frameRateNew == s.frameRate || isempty(s.frameRateNew)
        out = data;
    else
        out = resample(data, p, q);
    end
end
%%
function pupilFolder = parseFileName(filename, prefix)
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
    pupilFolder = [prefix y m2 d];
end
    


        
    


