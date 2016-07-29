function makeImageUserMenu
	global state gh
	
	if ishandle(state.internal.userSettingsMenu)
		delete(state.internal.userSettingsMenu);
	end
	
	if ~isempty(state.userSettingsPath)
		state.internal.userSettingsMenu=uimenu(gh.mainControls.figure1, 'Label', 'Users');
		flist=dir(fullfile(state.userSettingsPath, '*.usr'));
		uimenu(state.internal.userSettingsMenu, 'Label', state.userSettingsPath, 'Enable', 'on');
		
		for counter=1:length(flist)	
			if counter==1
				uimenu(state.internal.userSettingsMenu, 'Label', flist(counter).name, 'Callback', 'selectImageUserFromMenu' ...
					, 'Separator', 'on');
			else
				uimenu(state.internal.userSettingsMenu, 'Label', flist(counter).name, 'Callback', 'selectImageUserFromMenu');
			end
		end
	end		
	
	