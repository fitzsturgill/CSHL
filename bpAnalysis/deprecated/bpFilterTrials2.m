function ti = bpFilterTrials2(TE, varargin)
% 5/4/16
    counter = 1;
    ti = ones(size((TE.TrialTypes))); % assuming that you'll always have TrialTypes
    while counter+1 <= length(varargin) 
        field = varargin{counter};
        val = varargin{counter+1};
        theseMatches = ismember(TE.(field), val);
        if isempty(theseMatches)
            disp('wtf');
        end
        ti = ti & theseMatches;
        counter=counter+2;
    end