function rmsNoise = CSC_RMS(varargin)
    defaults = {...
        'filepath', [];... % if empty, interactive selection
        'channels', 1:32;...
        'indexRange', [1 1875];...
        'timestampRange', [];... % if time stamp range specified, overrides indexRange
        'Fs', 32000;...
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings

    if isempty(s.filepath)
        [fname, pname] = uiputfile('path', 'Choose CSC path...');
        if pname == 0
            rmsNoise = [];
            return
        else
            filepath = pname;
            disp(pname);
        end
    else
        filepath = s.filepath;
    end
       
    rmsNoise = NaN(numel(s.channels), 1);
    h = waitbar(0, 'RMS noise');   
    for channel = s.channels
        cfilename = sprintf('CSC%d.ncs', channel);
        try
            if isempty(s.timestampRange) && ~isempty(s.indexRange)
                [Samples, header] = Nlx2MatCSC(fullfile(filepath, cfilename),[0 0 0 0 1],1,2,s.indexRange);
            elseif ~isempty(s.timestampRange)
                [Samples, header] = Nlx2MatCSC(fullfile(filepath, cfilename),[0 0 0 0 1],1,4,s.timestampRange);
            else
                [Samples, header] = Nlx2MatCSC(fullfile(filepath, cfilename),[0 0 0 0 1],1,1);
            end
            Samples = reshape(Samples, numel(Samples), 1);
            ADBitVolts = sscanf(header{16}, '-ADBitVolts %f');
            Samples = Samples * ADBitVolts * 1e6;
            
            passBand = [10 10000]; % Hz
            [b, a] = butter(5, passBand/s.Fs * 2);
            filtData = filtfilt(b,a,Samples); 
            rmsNoise(channel) = rms(filtData);
        catch
        end
        waitbar(channel/32);
    end
    close(h);
end
