function [varargout] = loadAverage(filename, wavename)

if nargin == 0
	[fname, pname] = uigetfile('*avg.mat', 'Choose Save Path and Name for Wave object...');
	if isnumeric(pname)
		return
	else
		filename = [pname fname];
	end
end    

[path,fname,ext]=fileparts(filename);

if isempty(ext)
    ext = '.mat';
end

if isempty(path)
    path=pwd;
end

out=load([path '\' fname ext],'-mat');
name=char(fieldnames(out));
out=getfield(out,name);

if nargin == 2
    name=wavename;
end

if iswave(name)
    error(['loadWave: wave ' name ' already exists.  Use loadWaveo to overwrite']);
end

waveo(name, out.data, ...
		'xscale', out.xscale, ...
		'yscale', out.yscale, ...
        'zscale', out.zscale,...
		'UserData', out.UserData, ...
		'note', out.note, ...
		'timeStamp', out.timeStamp);

if nargout == 1
	varargout{1}=out;
end


%% now load components

components = out.UserData.Components;

for i = 1:length(components)
    filename = [components{1, i} '.mat'];
    loadWaveo(filename);
end