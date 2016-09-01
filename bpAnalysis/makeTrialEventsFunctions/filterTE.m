function validTrials = filterTE(TE, varargin) 
    
%     validTrials = zeros(length(TE.trialNumber), 1);
    % parse input parameter pairs
    counter = 1;
    while counter+1 <= length(varargin) 
        prop = varargin{counter};
        val = varargin{counter+1};
        
        if counter == 1
            validTrials(:,1) = ismember(TE.(prop), val); % initialize here
        else
            thisComparison = ismember(TE.(prop), val);
            if size(thisComparison, 2) > 1
                thisComparison = thisComparison';
            end
            validTrials = validTrials & thisComparison;
        end        
        counter=counter+2;
    end