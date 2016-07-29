function timerCallPackageFunctions(type, package)
	global state gh
	
	if nargin==2
		if iscell(package)
			pList=package;
		elseif isnumeric(package) & ~isempty(package)
			pList=state.timer.packageList(package);
		elseif ischar(package)
			pList=state.timer.packageList(strcmp(state.timer.packageList, package));
		else
			disp('timerCallPackageFunctions: package is of unknown type');
		end
	else
		pList=state.timer.packageList;
	end
	
	if strcmp(type, 'Abort')
		state.timer.abort=1;
		pList=state.timer.packageList;
	end
	
	for counter=1:length(pList)
		if state.timer.activePackages(timerPackageIndex(pList{counter})) & exist(['timer' type '_' pList{counter} '.m'])==2
		%	disp(['CALLING:  timer' type '_' pList{counter} ';']);
			eval(['timer' type '_' pList{counter} ';']);
		%	disp(['           Done']);
		end
	end
	
	if strcmp(type, 'Abort')
		timerCheckIfAllAborted
	end