function coh = bpCalcCoherence(data1, data2, Fs, varargin)
    %% compute coherence between 2 signals 
    % data 1 and 2-   (in form samples x trials)
    %% optional parameters, first set defaults
    defaults = {...
        'trialave', false;... 
        'err', [2 0.05];...
        'tapers', [3 5];... 
        'pad', true;... 
        'fpass', [0 20];...
        'whiten', false;...    % 1 = calculate cross spectrogram of temporal derivatives of signals, other options reserved for future (e.g. detrend, zscore, etc.)
        'uniformOutput', true;...            % not currently implemented, idea is to set to 0 if acqs are going to be variable in length (data in cell array)
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    if ~all(size(data1) == size(data2))
        error('data must be equal in size');
    end
    
    switch s.whiten
        case true
            data1 = diff(data1, 1, 1);
            data2 = diff(data2, 1, 1);
        otherwise
    end
    
    params.Fs = Fs;
    params.trialave = s.trialave;
    params.err = s.err;
    params.tapers = s.tapers;
    params.pad = s.pad;
    params.fpass = s.fpass;

    coh = struct(...
        'C', [],...
        'phi', [],...
        'S12', [],...
        'S1', [],...
        'S2', [],...    
        'f', [],...
        'confC', [],...
        'phistd', [],...
        'Cerr', [],...
        'settings', s,...
        'Fs', Fs...
        );

    [coh.C, coh.phi, coh.S12, coh.S1, coh.S2, coh.f, coh.confC, coh.phistd, coh.Cerr] = coherencyc(data1, data2, params);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    