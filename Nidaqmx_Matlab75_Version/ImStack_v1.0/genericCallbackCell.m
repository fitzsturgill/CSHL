function genericCallbackCell(handle)
global gh state
	
	value = get(gh.imageProcessingGUI.fileName, 'Value');
	genericCallback(handle);
	fieldname = getUserDataField(handle, 'Global');
	
	[structureName,structureToEnd,lastField] = structNameParts(fieldname);
	eval([structureToEnd '.cell.' lastField '{value} = ' fieldname ';']);
	
	