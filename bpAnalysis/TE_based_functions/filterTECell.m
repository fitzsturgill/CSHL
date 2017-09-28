function validTrials = filterTECell(TE, varargin) 
    
%     validTrials = zeros(length(TE.trialNumber), 1);
    % parse input parameter pairs
    cellPos = 1;  % in future make this flexible so you can for instance test if last event or state is such and such
    counter = 1;
    while counter+1 <= length(varargin) 
        prop = varargin{counter};
        val = varargin{counter+1};
        
        if counter == 1
            validTrials = cellfun(@(x) ismember(x(1), val), TE.(prop)); % initialize here
            if isnan(val) ==1
                validTrials = cellfun(@(x) isnan(x(1)), TE.(prop));
            end
        else
            validTrials = validTrials & cellfun(@(x) ismember(x(1), val), TE.(prop));
             if isnan(val) ==1
                validTrials = validTrials & cellfun(@(x) isnan(x(1)), TE.(prop));
            end
        end
        
        counter=counter+2;
    end