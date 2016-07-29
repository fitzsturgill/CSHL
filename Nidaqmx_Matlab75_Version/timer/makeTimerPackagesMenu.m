function makeTimerPackagesMenu
	global state gh
	
	children=get(gh.timerMainControls.Packages, 'Children');
	if ~isempty(children)
		delete(children);
	end
	
	if ~isempty(state.timer.packagesPath)
		flist=dir(fullfile(state.timer.packagesPath, 'timerInit_*.m'));
		uimenu(gh.timerMainControls.Packages, 'Label', state.timer.packagesPath, 'Enable', 'off', 'Callback', 'setPackagesPath');
		
		for counter=1:length(flist)	
			if counter==1
				uimenu(gh.timerMainControls.Packages, 'Label', flist(counter).name(11:end-2), 'Callback', 'selectPackageFromMenu' ...
					, 'Separator', 'on');
				state.timer.packageList={flist(counter).name(11:end-2)};
			else
				uimenu(gh.timerMainControls.Packages, 'Label', flist(counter).name(11:end-2), 'Callback', 'selectPackageFromMenu');
				state.timer.packageList=[state.timer.packageList {flist(counter).name(11:end-2)}];
			end
		end
			
		state.timer.activePackages=zeros(1,length(flist));
		state.timer.packageStatus=zeros(1,length(flist));
		state.timer.initializedPackages=zeros(1,length(flist));
		state.timer.pausedPackages=zeros(1,length(flist));
    else
        disp(['makeTimerPackagesMenu : Error: Packages path ' state.timer.packagesPath ...
            ' not found']);
		uimenu(gh.timerMainControls.Packages, 'Label', 'Set package path...', 'Enable', 'on', 'Callback', 'setPackagesPath');
		state.timer.activePackages=[];
		state.timer.initializedPackages=[];
		state.timer.packageStatus=[];
		state.timer.packageList={};
		state.timer.pausedPackages=[];
	end		