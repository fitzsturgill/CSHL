function st = b_backup( do_all )

% Fitz notes-   should make user-specific preferences for work_path and
% dest_path directories.  I should detect if dest_path directory exists
% rather than the weird success flag I use to determine whether an update
% has been succesful below.  I should also store backup file (time stamp of
% last backup) in the destination directory rather than the source
% directory (as opposed to the preferences which can be stored locally in
% the matlab path

% alex vaughn authorship?
%BACKUP   Creates backup files for Matlab work.
%   STATUS = BACKUP(LAST_BACKUP) saves all files in work directory that
%   have been modified since last backup to a destination directory
%   specified in the program code. ST = 1 if all files were successfully
%   copied and 0 otherwise.
%
%   See also FINISH.

% Creating machine independent work path
% work_path = fullfile(matlabroot,'work\Balazs\');
% work_path = 'C:\Bpod\Data';
% work_path = 'C:\Users\Adam\Documents\Repos\Bpod_Fitz\Data';
work_path = 'C:\Users\Adam\BpodUser\Data';

% Destination
% global DATAPATH
% destroot = [DATAPATH 'Matlab work backup\Balazs\'];
% destroot = '\\science\Kepecs\Balazs\Matlab_backup\';
% destroot = 'E:\FitzRig2\Data';
destroot = 'Z:\FitzRig2\Data';
if nargin == 0,
    do_all = 0;
end

% If doing incremental update, load last_backup timestamp
if do_all == 0
    try
        cd(work_path)
        load last_backup
    catch e
        if strcmp(e.identifier,'MATLAB:load:couldNotReadFile')
            fprintf('ERROR :: Could not find last_backup file.  NOT BACKING UP')
            return
        end
    end
    
    fprintf('Backing up all files since %s\n',datestr(last_backup,31))
    
else
    last_backup = 0;
end

    
% Copy
st = 1;   % status variable
dir_list{1} = work_path;
succeeded = false;
while ~isempty(dir_list)
    files = dir(dir_list{1});
    files = files(3:end);
    lf = length(files);
    for i = 1:lf
        if ~files(i).isdir
            % Not a directory
            dv = datevec(files(i).date);
            df = dv - last_backup;
            df(end+1) = -1;
            fdf = find(df);
            if df(fdf(1)) > 0
                source = fullfile(dir_list{1},files(i).name);
                destspec = dir_list{1}(length(work_path)+1:end);
                dest = fullfile(destroot,destspec);
                cmnd = ['copyfile(''',source,''',''',dest,''');'];
                fprintf('Copying %s\n     to %s \n',source,dest)

                try
                    eval(cmnd);
                    succeeded = true;
                catch
                    st = 0;
                    disp(['Unable to copy: ' source]);
                    lasterror
                end
            end
        else
            % A directory
            %put new directory on the directory list
            new_dir = fullfile(dir_list{1},files(i).name);
            dirspec = dir_list{1}(length(work_path)+1:end);
            ff = fullfile(destroot,dirspec);
            cmnd = ['status = mkdir(''',ff,''',''',files(i).name,''');'];
            eval(cmnd);
            if isempty(find(strcmp(dir_list,new_dir)))
                dir_list{end+1} = new_dir;
            end
        end
    end
    dir_list(1) = [];
end
%%
if succeeded
    cd(work_path)
    last_backup = clock
    save last_backup last_backup
else
    disp('BACKUP FAILED');
end