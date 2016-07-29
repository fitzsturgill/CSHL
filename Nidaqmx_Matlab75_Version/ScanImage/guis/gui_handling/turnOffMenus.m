function turnOffMenus

	global gh state
	turnoffPullDownMenu(gh.standardModeGUI.Settings, 'Edit Configuration...');
	turnoffPullDownMenu(gh.standardModeGUI.Settings, 'Channels...');
	turnoffPullDownMenu(gh.standardModeGUI.Settings, 'Get PMT Offsets...');
	turnoffPullDownMenu(gh.standardModeGUI.File, 'Load User Settings...');
	turnoffPullDownMenu(gh.standardModeGUI.File, 'Load Configuration...');

	if ishandle(state.internal.userSettingsMenu)
		 set(state.internal.userSettingsMenu, 'Enable', 'off');
	end
	if ishandle(state.internal.configurationsMenu)
		 set(state.internal.configurationsMenu, 'Enable', 'off');
	end
%	set(get(gh.standardModeGUI.figure1, 'children'), 'Enable', 'Off');
