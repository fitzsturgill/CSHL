function cxcg = bpCalcSpectrogram(data1, Fs, varargin)
    %% compute  spectrogram 
    % -   (in form samples x trials)
    %% optional parameters, first set defaults
    defaults = {...
        'trialave', false;... % 8/28/2016- changed channels default from [] to 1
        'err', [2 0.05];...
        'tapers', [5 9];... 
        'pad', true;... 
        'fpass', [0 20];...
        'movingwin', [5 1];...
        'uniformOutput', 1;...            % not currently implemented, idea is to set to 0 if acqs are going to be variable in length (data in cell array)
        'whiten', false;...      % 1 = calculate cross spectrogram of temporal derivatives of signals, other options reserved for future (e.g. detrend, zscore, etc.)
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings

    switch s.whiten
        case 1
            data1 = diff(data1, 1, 1);
        otherwise
    end


    params.Fs = Fs;
    params.trialave = s.trialave;
    params.err = s.err;
    params.tapers = s.tapers;
    params.pad = s.pad;
    params.fpass = s.fpass;

    cxcg = struct(...
        'S', [],...
        't', [],...    
        'f', [],...
        'settings', s,...
        'Fs', Fs...
        );


    [cxcg.S, cxcg.t, cxcg.f] = mtspecgramc(data1, s.movingwin, params);