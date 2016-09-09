function TE = addPupilometryToTE(TE, rootPath, varargin)

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
        'rootPath', [];...
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
    nFrames = round(s.frameRate * s.duration); % should be an integer anyway
    matFields = {'eyeArea', 'blinkDetected',...
        'pupArea', 'pupDiameter', 'pupResidual'};
    pupil = struct();
    pupil.loadSettings = s; % scalar
    pupil.settings = cell(length(TE.filename), 1);
    pupil.xData = 0:dX:(nFrames - 1) * dX;
    for field = matFields
        pupil.(field{:}) = NaN(length(TE.filename), nFrames); % fill with NaNs
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
            pupil.settings{tei} = loaded.pupilData.settings;          
            pupil.eyeArea(tei, 1:framesToLoad) = loaded.pupilData.eye.area(1:framesToLoad);
            pupil.blinkDetected(tei, 1:framesToLoad) = loaded.pupilData.eye.blinkDetected(1:framesToLoad);
            pupil.pupArea(tei, 1:framesToLoad) = loaded.pupilData.pupil.area(1:framesToLoad);
            pupil.pupDiameter(tei, 1:framesToLoad) = loaded.pupilData.pupil.diameter(1:framesToLoad);
            pupil.pupResidual(tei, 1:framesToLoad) = loaded.pupilData.pupil.circResidual(1:framesToLoad);            
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
    
    
        
    


