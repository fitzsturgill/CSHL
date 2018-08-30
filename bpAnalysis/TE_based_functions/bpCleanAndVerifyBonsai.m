function bpCleanAndVerifyBonsai(readonly)
    if nargin < 1
        readonly = 0;
    end
%     Function to delete weird spurious files that Bonsai spit out, there are
%     extra files for certain triggers that I can detect because their
%     datemodified timestamp is < n seconds from the previous timestamp
%     which is impossible given that Bpod trials are at least 11 seconds
%     long + the ITI

%     Also verifies that the number of Pupil.mat, Pupil.avi, Whisk.avi, and
%     Whisk.diff files matches the number of trials in the Bpod session
%     data
    shortestITI = 8; 
   
    sessions = bpLoadSessions; % load the session file
    if isempty(sessions)
        return
    end     
    for scounter = 1:length(sessions)
        fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n') 
        spath = sessions(scounter).filepath;
        sname = sessions(scounter).filename;
        nTrials = sessions(scounter).SessionData.nTrials;
        teDelta = diff(sessions(scounter).SessionData.TrialStartTimestamp);
        fileTypes = {'Pupil_', '*upil_*.avi';... 
                    'Whisk_', '*hisk_*.avi';...
                    'WhiskDiff_', '*hiskDiff_*.csv';...
                    };
        remainingFiles = zeros(size(fileTypes, 1), 1);
        remainingMatFiles = 0;
        combined = 0; % to be set to 1 in case bonsai files are in the Combined_ folder
        for counter = 1:size(fileTypes, 1)
            if ~combined
                targetFolder = parseFileName(sname, fileTypes{counter, 1});
                targetFolder = fullfile(spath, targetFolder);
                if ~exist(targetFolder, 'dir') % it's a folder    
                    combined = 1;
                end
            end

            if combined
                targetFolder = parseFileName(sname, 'Combined_');
                targetFolder = fullfile(spath, targetFolder);
                if ~exist(targetFolder, 'dir') % it's a folder    
                    disp('*** target folders are not found ***');
                    return                
                end            
            end

            cd(targetFolder);
            fs = dir(fileTypes{counter, 2});
            fileList = {};
            [fileList{1:length(fs)}] = fs(:).name; % see deal documentation I think...
            [fileList,ix] = sort_nat(fileList); % alphanumeric sorting
            dmDelta = seconds(diff(datetime({fs(ix).date})));

            tomove = fileList(find(dmDelta < shortestITI) + 1);
            nRemaining = sum(dmDelta >= shortestITI) + 1; % remember to add on + 1 for the first file...
            if ~isempty(tomove)
                backupFolder = fullfile(targetFolder, ['Spurious' filesep]);
                ensureDirectory(backupFolder);
                for fcounter = 1:length(tomove)
                    if ~readonly
                        movefile(fullfile(targetFolder, tomove{fcounter}), fullfile(backupFolder, tomove{fcounter}));
                        fprintf('*** Moved spurious %s to folder: %s\n', fullfile(targetFolder, tomove{fcounter}), backupFolder)
                    end
                    if counter == 1 % also try moving the Pupil_*.mat file if necessary
                        [~, matName, ~] = fileparts(tomove{fcounter});
                        if exist(fullfile(targetFolder, matName, '.mat'), 'file')
                            if ~readonly
                                movefile(fullfile(targetFolder, matName, '.mat'), fullfile(backupFolder, tomove{fcounter}));
                                fprintf('*** Also moved %s \n', fullfile(matName, '.mat'))
                            end
                        end

                    end               
                end
            end  
            % get remaining number of .mat files for verification
            if counter == 1
                fsMat = dir('*upil_*.avi');
                nMatFiles = length(fsMat);          
                fprintf('   Bonsai file report for %s :\n', sname)
                fprintf('   Trials = %d\n', nTrials)
                fprintf('   Pupil_*.mat = %d, difference of %d\n', nMatFiles, nMatFiles - nTrials);        
            end
            fprintf('   %s = %d, difference of %d\n', fileTypes{counter, 2}, nRemaining, nRemaining - nTrials);
        end
    end
    fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')    




end

function targetFolder = parseFileName(filename, prefix)
    % extract date and generate file name matching for example
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
    targetFolder = [prefix y m2 d filesep];
end
    
    
    