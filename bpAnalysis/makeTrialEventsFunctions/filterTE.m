function validTrials = filterTE(TE, varargin) 
    
%     validTrials = zeros(length(TE.trialNumber), 1);
    % parse input parameter pairs
    counter = 1;
    while counter+1 <= length(varargin) 
        prop = varargin{counter};
        val = varargin{counter+1};
        
        if counter == 1
            validTrials = ismember(TE.(prop), val); % initialize here
        else
            validTrials = validTrials & ismember(TE.(prop), val);
        end
        
        counter=counter+2;
    end