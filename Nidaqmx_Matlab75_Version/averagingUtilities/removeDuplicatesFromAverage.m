function removeDuplicatesFromAverage
    % FS 8/21/2013
	[fname, pname] = uigetfile('*.mat', 'Select Average');
	if isnumeric(pname)
		return
	else
		filename = [pname fname];
    end
    
    ext='';
    [path,fname,ext]=fileparts(filename);
    if isempty(ext)
        ext='.mat';
    end

    if isempty(path)
        path=pwd;
    end

    out=load([path '\' fname ext],'-mat');
    name=char(fieldnames(out));
    out=getfield(out,name);
    
    waveo(name, out.data, ...
            'xscale', out.xscale, ...
            'yscale', out.yscale, ...
            'zscale', out.zscale,...
            'UserData', out.UserData, ...
            'note', out.note, ...
            'timeStamp', out.timeStamp);    
    
    
    
    components = avgComponentList(name);
    
    uniqueComp = unique(components);
    
    
    setWaveUserDataField(name, 'nComponents', length(uniqueComp));
    setWaveUserDataField(name, 'Components', uniqueComp);

    saveWave(name, filename);
    kill(name);
    
    
    
    
    