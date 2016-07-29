function timerSetActiveStatus(package, status)
	if nargin<2
		status=1;
	end
	if nargin<1
		error('timerSetActiveStatus: please provide package name');
	end;
	index=timerPackageIndex(package);
	if isempty(index)
		error(['timerSetActiveStatus: unknown package: ' package '--> Try correcting package path in timer.ini']);
	end
	
	global state gh
	state.timer.activePackages(index)=status;
	if status	% package was turned on;  if not initialized, then do it now
		if ~state.timer.initializedPackages(index)
			timerCallPackageFunctions('Init', index);
			state.timer.initializedPackages(index)=1;
		end	
	end
	if ishandle(gh.timerMainControls.Packages)	% set the flag in the menu
		menuIndex=length(state.timer.packageList)-index+1;
		children=get(gh.timerMainControls.Packages, 'Children');
		if status
			set(children(menuIndex), 'Checked', 'on');
		else
			set(children(menuIndex), 'Checked', 'off');
		end
	end
	