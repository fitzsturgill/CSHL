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
    for counter = 1:length(sessionnames)
        sessionname = sessionnames{counter};
        pupilFolder = parseFileName(sessionname); % see subfunction
        pupilPath = fullfile(rootPath, [pupilFolder s.folderSuffix{counter}], filesep); % filesep returns system file separator character
        if ~isdir(pupilPath)
            warning(['*** pupil directory: ' pupilPath ' does not exist ***']);
            continue
        end
        cd(pupilPath);
        fs = dir('*upil_*.mat');
        if length(fs) == 0
            warning(['*** No Pupil_ files found in ' pupilFolder]); % you need to analyze it first!
            continue;
        end
        fileList = {};
        [fileList{1:length(fs)}] = fs(:).name; % see deal documentation I think...
        [fileList,~] = sort_nat(fileList); % alphanumeric sorting
        
        % find matching session indices within TE
        si = find(filterTE(TE, 'filename', sessionname));
        for i = 1:length(si)
            tei = si(i); % TE index number
            % make sure numbers match
            fname = fileList{i}; % there should exist a .mat file for every index in TE
            pupNumber = str2double(fname(strfind(fname, '_') + 1:strfind(fname, '.mat') - 1)); % Pupil_11.mat, extract 11
            if TE.trialNumber(tei) ~= pupNumber + 1; % bonsai currently names saved videos starting at 0 in my pipeline
                warning(['*** pupil file numbering mismatch for session ' sessionname ' skipping... ***']);
                break
            end
            loaded = load(fileList{i});
            framesAvailable = loaded.pupilData.currentFrame(end);
            framesToLoad = min(nFrames, max(loaded.pupilData.currentFrame)); % bonsai seems most often to add one extra frame, still handling case where a few frames are dropped, see below
            newFrames = round(framesToLoad * (s.frameRateNew/s.frameRate)); % is round correct? 
            pupil.frameRate(tei,1) = s.frameRateNew;
            pupil.settings{tei} = loaded.pupilData.settings;          
            pupil.eyeArea(tei, 1:newFrames) = changeRate([loaded.pupilData.eye.area(1:framesToLoad) repmat(loaded.pupilData.eye.area(end), 1, nFrames - framesToLoad)]);
            pupil.blinkDetected(tei, 1:newFrames) = changeRate([loaded.pupilData.eye.blinkDetected(1:framesToLoad) repmat(loaded.pupilData.eye.blinkDetected(end), 1, nFrames - framesToLoad)]);
            pupil.pupArea(tei, 1:newFrames) = changeRate([loaded.pupilData.pupil.area(1:framesToLoad) repmat(loaded.pupilData.pupil.area(end), 1, nFrames - framesToLoad)]);
            pupil.pupDiameter(tei, 1:newFrames) = changeRate([loaded.pupilData.pupil.diameter(1:framesToLoad) repmat(loaded.pupilData.pupil.diameter(end), 1, nFrames - framesToLoad)]);
            pupil.pupResidual(tei, 1:newFrames) = changeRate([loaded.pupilData.pupil.circResidual(1:framesToLoad) repmat(loaded.pupilData.pupil.circResidual(end), 1, nFrames - framesToLoad)]);
            pupil.startTime(tei, 1) = TE.(s.startField){tei}(1);

            if counter == 1 && i == 1
                startX = TE.(s.startField){1}(1) - TE.(s.zeroField){1}(1);
                pupil.xData = pupil.xData + startX;
                blStartP = 1; % just start at beginning of video for baseline
                blEndP = bpX2pnt(0, s.frameRateNew, startX); % go to zero point        
            end
        end  
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

    %% nested function
    function out = changeRate(data)
        if s.frameRateNew == s.frameRate || isempty(s.frameRateNew)
            out = data;
        else
            [p, q] = rat(s.frameRateNew/s.frameRate);
            out = resample(data, p, q);
        end
    end
        
end
%%
function pupilFolder = parseFileName(filename)
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
    pupilFolder = ['Pupil_' y m2 d];
end
    


        
    


