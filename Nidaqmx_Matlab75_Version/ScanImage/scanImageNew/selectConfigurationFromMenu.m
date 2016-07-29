function selectConfigurationFromMenu(forceName)
	global state gh

	h=gcbo;
		
	setStatusString('Reading config file...');
	children=get(gh.standardModeGUI.Configurations, 'Children');

	if nargin<1
		index=find(children==h);
	else
		index=find(strcmp(get(children, 'Label'), forceName));
	end
	[pathstr,name,ext,versn] = fileparts(get(children(index), 'Label'));
	state.configName=name;
	state.configPath=get(children(end), 'Label');

	loadStandardModeConfig;
