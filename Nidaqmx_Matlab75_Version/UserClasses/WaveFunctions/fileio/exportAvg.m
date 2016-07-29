function exportAvg(waveList, outputFile)

	if any(iswave(waveList))
		if ~iscell(waveList) & ~ischar(waveList) 
			waveList=inputname(1);
		end
	else
		error('exportWave : expect wave name or cell array of wave names as input');
	end
	
	if ischar(waveList)
		waveList={waveList};
	end

	if nargin==1
		[filename, pathname] = uiputfile('*.itx', 'Save wave as');
		if isempty(filename)
			return
		end
		outputFile=fullfile(pathname, filename);
	end
	if isempty(find(outputFile=='.'))
		outputFile=[outputFile '.itx'];
	end

	[fid, message] = fopen(outputFile, 'w');
	if fid==-1
		message
		error(['exportAvg : could not open ' outputFile ' for output']);
	end
	
	fprintf(fid, 'IGOR\n');
	for avgName=waveList
		disp([' Exporting : ' avgName{1} ' and its components']);
		for name=[avgName avgComponentList(avgName{1})];
			if iswave(name{1})
				fprintf(fid, 'WAVES %s\nBEGIN\n', name{1});
				data=get(name{1}, 'data');
				for counter=1:size(data, 2)
					fprintf(fid, '   %e\n', data(counter));
				end
				fprintf(fid, 'END\n');	
				xscale=get(name{1}, 'xscale');
				fprintf(fid, 'X SetScale/P x %f, %f, "", %s\n', xscale(1), xscale(2), name{1});
				if ~isempty(getWaveUserDataField(name{1}, 'nComponents'))
					cList=(getWaveUserDataField(name{1}, 'Components'));
					cText=['X Note ' name{1} ', "nAvg=' num2str(getWaveUserDataField(name{1}, 'nComponents')) ';parts='];
					for comp=cList
						cText=[cText comp{1} ','];
					end
					cText=[cText(1:end-1) ';"'];
					fprintf(fid, '%s\n\n', cText);
				end
				fprintf(fid, '\n');
			else
				disp(['exportAvg : ' name{1} ' does not exist or is not a wave.  Skipping...']);
			end
		end
	end
	fclose(fid);