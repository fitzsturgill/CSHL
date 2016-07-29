function turnOnMenus

	global gh state

	turnonPullDownMenu(gh.standardModeGUI.Settings, 'Edit Configuration...');
	turnonPullDownMenu(gh.standardModeGUI.Settings, 'Channels...');
	turnonPullDownMenu(gh.standardModeGUI.File, 'Load User Settings...');
	turnonPullDownMenu(gh.standardModeGUI.File, 'Load Configuration...');
	turnonPullDownMenu(gh.standardModeGUI.Settings, 'Get PMT Offsets...');
	if ishandle(state.internal.userSettingsMenu)
		 set(state.internal.userSettingsMenu, 'Enable', 'On');
	end
	if ishandle(state.internal.configurationsMenu)
		 set(state.internal.configurationsMenu, 'Enable', 'On');
	end	
	
%	set(get(gh.standardModeGUI.figure1, 'children'), 'Enable', 'On');
