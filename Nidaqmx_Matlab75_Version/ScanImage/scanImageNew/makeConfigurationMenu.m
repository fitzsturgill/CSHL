function makeConfigurationMenu
	global state gh
	
	children=get(gh.standardModeGUI.Configurations, 'Children');
	if ~isempty(children)
		delete(children);
	end
	
	if ~isempty(state.configPath)
		flist=dir(fullfile(state.configPath, '*.cfg'));
		uimenu(gh.standardModeGUI.Configurations, 'Label', state.configPath, 'Enable', 'on', 'Callback', 'setConfigurationsPath');
		
		for counter=1:length(flist)	
			if counter==1
				uimenu(gh.standardModeGUI.Configurations, 'Label', flist(counter).name, 'Callback', 'selectConfigurationFromMenu' ...
					, 'Separator', 'on');
			else
				uimenu(gh.standardModeGUI.Configurations, 'Label', flist(counter).name, 'Callback', 'selectConfigurationFromMenu');
			end
		end
	else		
		uimenu(gh.standardModeGUI.Configurations, 'Label', 'Set Configurations Path...', 'Enable', 'on', 'Callback', 'setConfigurationsPath');
	end	
	