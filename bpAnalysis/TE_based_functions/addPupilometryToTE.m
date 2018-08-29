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

    wtf.files = [];
    wtf.trials = [];
    wtf = repmat(wtf, length(sessionnames), 1);
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
%         fs = dir('*upil_*.mat');
        fs = dir('*upil_*.avi'); % use videos to get time stamps, then generate .mat file name later on which will be derived from the timstamped video file
        if length(fs) == 0
            warning(['*** No Pupil_ files found in ' pupilFolder]); % you need to analyze it first!
            continue
        end
        fileList = {fs(:).name};        
        [fileList,ix] = sort_nat(fileList); % alphanumeric sorting
        dmDelta = seconds(diff(datetime({fs(ix).date})));
        dmDelta = dmDelta(:);
        if any(dmDelta < 2)
            error('there are spurious pupil files that you need to delete or deal with');
        end
%         wtf(counter).files = dmDelta;
        % find matching session indices within TE either by matching numeric suffix or total number of files 
        si = find(filterTE(TE, 'filename', sessionname));
        teDelta = diff(TE.TrialStartTimestamp(si));
        teDelta = teDelta(:);
%         wtf(counter).trials = teDelta;
        %% verify correct numbering of pupil files
%         methods = [0 0]; % 2 methods of checking
        % method #1- see if first file is Pupil_0.mat file.....
        fname = fileList{1}; % there should exist a .mat file for every index in TE
        [~,fname,~] = fileparts(fname);
        firstNumber = str2double(fname(strfind(fname, '_') + 1:end));
%         if firstNumber == 0
%             methods(1) = 1;
%         end
        
        % note 180826 FS
        % sometimes there is 1 fewer bonsai acq, I think this is because
        % bonsai gets behind and I need to wait to stop bonsai after a
        % trial
        
        % note 180827 FS, sometimes there are way more bonsai acqs this is
        % because there are spurious triggers (is my pulse to bonsai too
        % long?) and you get nearly duplicate bonsai files saved 
%         if abs((length(si) - length(fileList))) > 1
%             methods(2) = 1;
%         end
        
        % correlate ITIs with date modified on pupil video files, try to
        % shift if necessary, assuming an alignment problem at the
        % beginning
        % cirshift on teDelta, if - works, then the first n trials don't have a pupil file (you need
        % to skip the first n trials)
        % otherwise, the first n pupil files don't have a trial ( you need
        % to skip the first n pupil files)
        
        goodShift = [];        
        ml = min(length(dmDelta), length(teDelta));
        maxRsq = 0;
        for shift = -3:3
            rsq = corr(dmDelta(1:ml), circshift(teDelta(1:ml), shift));
            maxRsq = max(rsq, maxRsq);
            if rsq > 0.9
                goodShift = shift + 1;
                fprintf('*** Rsq of trial and file ITIs = %.2f with shift of %d for %s ***\n', rsq, goodShift, sessionname); 
                break
            end
        end     

                                                                 
        if isempty(goodShift)  
            warning('*** pupil file alignment issue with max Rsq = %.2f for session %s skipping... *** \n', maxRsq, sessionname);
%             fprintf('*** First pupil suffix is %d, %d trials vs %d Pupil files for session %s ***, skipping\n', firstNumber, length(si), length(fileList), sessionname);            
            continue
        end
        
        %% loop through pupil .mat files from a given session
        wb = waitbar(0, num2str(counter));
        for i = 1:length(si) 
            tei = si(i); % TE index number    
            fi = i + goodShift;
            if fi <= 0 % there isn't a pupil file for this one
                continue
            end
           
            try
                [~,fname,~] = fileparts(fileList{fi});                 
                loaded = load([fname '.mat']);
            catch
                fprintf('*** %s didn''t exist in %s ***\n', [fname '.mat'], pupilPath)
                continue;
            end
            
            %% padding with end value
%             framesToLoad = min(nFrames, loaded.pupilData.currentFrame(end)); % bonsai seems most often to add one extra frame, still handling case where a few frames are dropped, see below
%             newFrames = s.duration * s.frameRateNew;
%             pupil.frameRate(tei,1) = s.frameRateNew;
%             pupil.settings{tei} = loaded.pupilData.settings;          
% 
%             pupil.eyeArea(tei, 1:newFrames) = changeRate([loaded.pupilData.eye.area(1:framesToLoad)' repmat(loaded.pupilData.eye.area(end), 1, nFrames - framesToLoad)]);
%             pupil.blinkDetected(tei, 1:newFrames) = changeRate([loaded.pupilData.eye.blinkDetected(1:framesToLoad)' repmat(loaded.pupilData.eye.blinkDetected(end), 1, nFrames - framesToLoad)]);
%             pupil.pupArea(tei, 1:newFrames) = changeRate([loaded.pupilData.pupil.area(1:framesToLoad)' repmat(loaded.pupilData.pupil.area(end), 1, nFrames - framesToLoad)]);
%             pupil.pupDiameter(tei, 1:newFrames) = changeRate([loaded.pupilData.pupil.diameter(1:framesToLoad)' repmat(loaded.pupilData.pupil.diameter(end), 1, nFrames - framesToLoad)]);
%             pupil.pupResidual(tei, 1:newFrames) = changeRate([loaded.pupilData.pupil.circResidual(1:framesToLoad)' repmat(loaded.pupilData.pupil.circResidual(end), 1, nFrames - framesToLoad)]);
%             pupil.startTime(tei, 1) = TE.(s.startField){tei}(1);
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

save(fullfile(rootPath, 'WTF.mat'), 'wtf');
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
    


        
    


