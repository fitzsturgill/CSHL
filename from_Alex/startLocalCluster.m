function startLocalCluster(force_restart)
%
% Starts a local parpool/matlabpool cluster of aporpriate size given teh
% number of cores on the current computer.  Compatible with OSX and
% WINDOWS, as well as matlabpool/parpool functions.
%
%
% Alex Vaughan, 2015

if nargin == 0
    force_restart = 0;
end

% Define max pool size for each architecture
switch computer
    case 'GLNXA64' % Unix!
        [~,max_pool_size] = unix('nproc');
        max_pool_size = max_pool_size - 2;
    case 'MACI64'  % OSX!
        [~,max_pool_size] = unix('sysctl -n hw.ncpu'); % Return value is char.
        max_pool_size = str2num(max_pool_size) - 2;
    otherwise
        error('Architecture not supported.')
end

% Start cluster (but don't restart an existing cluster)
if exist('parpool','file'),
    % For ~2014 and above
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        poolsize = 0;
    else
        poolsize = poolobj.NumWorkers;
    end
    if poolsize ~= max_pool_size || force_restart
        fprintf('Setting up pool size via parpool: %g.\n',max_pool_size)
        try
            delete(gcp('nocreate'));
        catch
            disp('Failed to delete parpool object')
        end
        parpool('local',max_pool_size);
    else
        fprintf('Pool size via parpool: %g.\n',max_pool_size)
    end
else
    % For ~2013 and below.
    fprintf('Setting up pool size via matlabpool: %g.\n',max_pool_size)
    if matlabpool('size') ~= max_pool_size || force_restart %#ok<*DPOOL>
        try
            matlabpool close
        end
        matlabpool max_pool_size;
    end
end