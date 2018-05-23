function whisk = addWhiskingToTE(TE, varargin)

% Adds whisking data to TE.  Pupil data for every session loaded into TE
% assumed to to be saved in a folder (e.g. \Whisking_160814\). Whisking data .csv
% files (e.g. Whisking_0.csv) are zero indexed
 
% TE: Trial Event structure containing fields trialNumber and filename
% rootPath: path containing session data and Whisking session folders

    %% optional parameters, first set defaults
    defaults = {...
        'sampleRate', 60;... % integer
        'sampleRateNew', 20;... 
        'duration', 11;...  % assumes constant trial duration across TE
        'zeroField', 'Us';... 
        'startField', 'PreCsRecording';... % extract time stamp w.r.t. start of Bpod trial
        'rootPath', [];...
        'normMode', 'bySession';...  % [bySession, byTrial]
        'folderPrefix', 'Combined_';...
        'folderSuffix', '';...
        'filePrefix', 'WhiskDiff_';...
        };
    
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings    
    if isempty(s.sampleRateNew)
        s.sampleRateNew = s.sampleRate;
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
    
    
    %% initialize pupil structure
    dX = 1/s.sampleRateNew;
    startX = []; % to be determined
    nSamples = s.sampleRate * s.duration; 
    normFields = {...
        'whisk', 'whiskNorm';...
        };
    whisk = struct();
    whisk.loadSettings = s; % scalar


    whisk.xData = linspace(0, s.duration - dX, s.duration * s.sampleRateNew);
    nsamplesNew = length(whisk.xData);
    whisk.startTime = NaN(length(TE.filename), 1);
    whisk.sampleRate = NaN(length(TE.filename), 1);
    
    for counter = 1:numel(normFields)
        whisk.(normFields{counter}) = NaN(length(TE.filename), nsamplesNew); % fill with NaNs
    end
    
    
    %% extract session dates in TE and generate matching whisk folder names
    sessionnames = unique(TE.filename);
    if ischar(s.folderSuffix)
        s.folderSuffix = repmat({s.folderSuffix}, length(sessionnames));
    end
    for counter = 1:length(sessionnames)
        sessionname = sessionnames{counter};
        whiskFolder = parseFileName(sessionname, s.folderPrefix); % see subfunction
        whiskPath = fullfile(rootPath, [whiskFolder s.folderSuffix{counter}], filesep); % filesep returns system file separator character
        if ~isdir(whiskPath)
            warning(['*** whisk directory: ' whiskPath ' does not exist ***']);
            continue
        end
        cd(whiskPath);
        fs = dir([s.filePrefix '*.csv']);
        if length(fs) == 0
            warning(['*** No ' s.filePrefix ' files found in ' whiskFolder]); % you need to analyze it first!
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
            whiskNumber = str2double(fname(strfind(fname, '_') + 1:strfind(fname, '.csv') - 1)); % whisk_11.mat, extract 11
            if TE.trialNumber(tei) ~= whiskNumber + 1 % bonsai currently names saved videos starting at 0 in my pipeline
                warning(['*** whisk file numbering mismatch for session ' sessionname ' skipping... ***']);
                break
            end
            % extract only first column using csv read
            loaded = dlmread(fileList{i}, ' ');
            loaded = loaded(:,1);
            samplesAvailable = size(loaded, 1);
            oldRate = samplesAvailable / s.duration;
            [p, q] = rat(s.sampleRateNew/oldRate);
            resampled = resample(loaded, p, q);
            whisk.sampleRate(tei, 1) = s.sampleRateNew;
            whisk.whisk(tei,:) = resampled';
            
            if counter == 1 && i == 1
                startX = TE.(s.startField){1}(1) - TE.(s.zeroField){1}(1);
                whisk.xData = whisk.xData + startX;
                blStartP = 1; % just start at beginning of video for baseline
                blEndP = bpX2pnt(0, s.sampleRateNew, startX); % go to zero point        
            end            
        end
            %% baseline subtract
        for i = 1:size(normFields, 1)
            rawField = normFields{i, 1};
            normField = normFields{i, 2};
            rawData = whisk.(rawField)(si, :);
            switch s.normMode
                case 'byTrial'
                    whisk.(normField)(si, :) = bsxfun(@rdivide, rawData, nanmean(rawData(:, blStartP:blEndP), 2));
                case 'bySession'
                    whisk.(normField)(si, :) = rawData / nanmean(nanmean(rawData(:, blStartP:blEndP), 2), 1); % divide by scalar
            end
        end
    end
end

function pupilFolder = parseFileName(filename, folderPrefix)
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
    pupilFolder = [folderPrefix y m2 d];
end