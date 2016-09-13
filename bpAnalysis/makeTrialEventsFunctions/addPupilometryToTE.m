function TE = addPupilometryToTE(TE, varargin)

% Adds pupil data to TE.  Pupil data for every session loaded into TE
% assumed to to be saved in a folder (e.g. \Pupil_160814\). Pupil data .mat
% files (e.g. Pupil_0.mat) are zero indexed

% TE: Trial Event structure containing fields trialNumber and filename
% rootPath: path containing session data and pupil session folders

    %% optional parameters, first set defaults
    defaults = {...
        'frameRate', 60;... % 8/28/2016- changed channels default from [] to 1
        'duration', 10;...  % assumes constant trial duration across TE
        'zeroField', 'Us';...
        'startField', 'PreCsRecording';... % extract time stamp w.r.t. start of Bpod trial
        'rootPath', [];...
        'normMode', 'bySession';...  % [bySession, byTrial]
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings    
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
    
    %% initialize pupil structure
    dX = 1/s.frameRate;
    startX = []; % to be determined
    nFrames = round(s.frameRate * s.duration); % should be an integer anyway
    normFields = {...
        'eyeArea', 'eyeAreaNorm';...
        'pupArea', 'pupAreaNorm';...
        'pupDiameter', 'pupDiameterNorm';...
        'pupResidual', 'pupResidualNorm';...
        };
    pupil = struct();
    pupil.loadSettings = s; % scalar
    pupil.settings = cell(length(TE.filename), 1);
    pupil.xData = 0:dX:(nFrames - 1) * dX;
    pupil.startTime = NaN(length(TE.filename), 1);
    pupil.frameRate = NaN(length(TE.filename), 1);
    pupil.blinkDetected = NaN(length(TE.filename), nFrames);
    
    for counter = 1:numel(normFields)
        pupil.(normFields{counter}) = NaN(length(TE.filename), nFrames); % fill with NaNs
    end
    
    %% extract session dates in TE and generate matching pupil folder names
    sessionnames = unique(TE.filename);

    for counter = 1:length(sessionnames)
        sessionname = sessionnames{counter};
        pupilFolder = parseFileName(sessionname); % see subfunction
        pupilPath = fullfile(rootPath, pupilFolder, filesep); % filesep returns system file separator character
        if ~isdir(pupilPath)
            warning(['*** pupil director: ' pupilPath ' does not exist ***']);
            continue
        end
        cd(pupilPath);
        fs = dir('Pupil_*.mat');
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
            framesToLoad = min(nFrames, max(loaded.pupilData.currentFrame));
            pupil.frameRate(tei,1) = s.frameRate;
            pupil.settings{tei} = loaded.pupilData.settings;          
            pupil.eyeArea(tei, 1:framesToLoad) = loaded.pupilData.eye.area(1:framesToLoad);
            pupil.blinkDetected(tei, 1:framesToLoad) = loaded.pupilData.eye.blinkDetected(1:framesToLoad);
            pupil.pupArea(tei, 1:framesToLoad) = loaded.pupilData.pupil.area(1:framesToLoad);
            pupil.pupDiameter(tei, 1:framesToLoad) = loaded.pupilData.pupil.diameter(1:framesToLoad);
            pupil.pupResidual(tei, 1:framesToLoad) = loaded.pupilData.pupil.circResidual(1:framesToLoad);
            pupil.startTime(tei, 1) = TE.(s.startField){1}(1);
            if counter == 1 && i == 1
                startX = TE.(s.zeroField){1}(1) - TE.(s.startField){1}(1);
                pupil.xData = pupil.xData - startX;
                blStartP = 1; % just start at beginning of video for baseline
                blEndP = bpX2pnt(0, s.frameRate, startX); % go to zero point        
            end
        end  
        % baseline subtract
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
    
    
        
    


