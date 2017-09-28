function validTrials = filterTE(TE, varargin) 
    % validTrials = nTrials x 1 logical array
    validTrials = ones(length(TE.filename), 1);
    testTrials = zeros(length(TE.filename), 1);
    % parse input parameter pairs
    counter = 1;
    while counter+1 <= length(varargin) 
        prop = varargin{counter};
        val = varargin{counter+1};
        testTrials(:,1) = ismember(TE.(prop), val);
        validTrials = validTrials & testTrials;
        counter=counter+2;
    end