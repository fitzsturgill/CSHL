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
    
    
    %% initialize whisk structure
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
    h = waitbar(0, 'Processing Whisking'); 
    
    %% assuming consistent trial durations and alignment
    startX = TE.(s.startField){1}(1) - TE.(s.zeroField){1}(1);
    whisk.xData = whisk.xData + startX;
    blStartP = 1; % just start at beginning of video for baseline
    blEndP = bpX2pnt(0, s.sampleRateNew, startX); % go to zero point        
    [p, q] = rat(s.sampleRateNew/s.sampleRate);
    for counter = 1:length(sessionnames)
        sessionname = sessionnames{counter};
        % try Whisk_ and Combined_ as prefixes
        whiskFolder = {parseFileName(sessionname, 'WhiskDiff_'), parseFileName(sessionname, 'Combined_')}; % see subfunction
        folderFound = 0;
        for j = 1:length(whiskFolder)
            whiskPath = fullfile(rootPath, [whiskFolder{j} s.folderSuffix{counter}], filesep); % filesep returns system file separator character
            if isdir(whiskPath)
                folderFound = 1; % whiskPath is correct....
                break
            end
        end
        if ~folderFound
            warning('*** Whisk directory does not exist ***');
            continue
        end
        cd(whiskPath);
        fs = dir([s.filePrefix '*.csv']);
        if length(fs) == 0
            warning(sprintf('*** No %s files found in %s or %s ***', s.filePrefix ,whiskFolder{1} ,whiskFolder{2})); 
            continue
        end
        fileList = {};
        [fileList{1:length(fs)}] = fs(:).name; % see deal documentation I think...
        [fileList,ix] = sort_nat(fileList); % alphanumeric sorting
        dmDelta = seconds(diff(datetime({fs(ix).date})));
        if any(dmDelta < 2)
            error('there are spurious whisk files that you need to delete or deal with');
        end        
        
        % find matching session indices within TE
        si = find(filterTE(TE, 'filename', sessionname));
        
        
        
        %% verify correct numbering of whisk files
        methods = [0 0]; % 2 methods of checking
        % method #1- see if first file is Pupil_0.mat file.....
        fname = fileList{1}; % there should exist a .mat file for every index in TE
        firstNumber = str2double(fname(strfind(fname, '_') + 1:strfind(fname, '.csv') - 1));
        if firstNumber == 0
            methods(1) = 1;
        end
        % note 180826 FS
        % sometimes there is 1 fewer bonsai acq, I think this is because
        % bonsai gets behind and I need to wait to stop bonsai after a
        % trial
        
        % note 180827 FS, sometimes there are way more bonsai acqs this is
        % because there are spurious triggers (is my pulse to bonsai too
        % long?) and you get nearly duplicate bonsai files saved         
        if abs((length(si) - length(fileList))) <= 1
            methods(2) = 1;
        end
        
        if ~any(methods)
%             warning(['*** pupil file numbering mismatch for session ' sessionname ' skipping... *** \n']);
            fprintf('*** First whisk suffix is %d, %d trials vs %d whisk files for session %s ***, skipping\n', firstNumber, length(si), length(fileList), sessionname);
            continue
        end
        
        for i = 1:min(length(si), length(fileList)) % in case we have one fewer pupil video
            tei = si(i); % TE index number
            % make sure numbers match
            fname = fileList{i}; % there should exist a .mat file for every index in TE

            % extract only first column using csv read
            try
                loaded = dlmread(fileList{i}, ' ');
            catch
                fprintf('*** Skipping %s from %s ***\n', fileList{i}, sessionname);
                continue
            end
            loaded = loaded(:,1);
            samplesAvailable = size(loaded, 1);            
            currentSamples = ceil(samplesAvailable * p/q);
            resampled = resample(loaded, p, q);
            whisk.sampleRate(tei, 1) = s.sampleRateNew;
            whisk.whisk(tei,1:min(nsamplesNew, currentSamples)) = resampled(1:min(nsamplesNew, currentSamples))';               
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
        waitbar(counter/length(sessionnames));
    end
    close(h);
end

function whiskFolder = parseFileName(filename, folderPrefix)
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
    whiskFolder = [folderPrefix y m2 d];
end