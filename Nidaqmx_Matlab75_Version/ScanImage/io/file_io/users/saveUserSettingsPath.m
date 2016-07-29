function saveUserSettingsPath(userPath)

	if nargin<1
		global state
		userPath=state.userSettingsPath;
	end
	save(fullfile(matlabroot, 'work', 'ScanImage', ['lastUserPath.mat']), ...
		'userPath', '-mat');
	

	