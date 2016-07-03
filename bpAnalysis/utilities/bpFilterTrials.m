function out = bpFilterTrials(SessionData, types, outcomes, filterField)
    % to not restrict based upon either types or outcomes, set to empty []
    % ex:   bpFilterTrials(SessionData, [1 2 4 5], []);   
    if nargin < 4
        filterField = 'IncludeTrials';
    end
    out = ones(1, SessionData.nTrials);
   
    % optionally include/exclude trials based upon filterField
    if isfield(SessionData, filterField)
        out = out & SessionData.(filterField);
    end
    
    if ~isempty(types)
        byType = ismember(SessionData.TrialTypes(1:SessionData.nTrials), types);
        out = out & byType;
    end
    if ~isempty(outcomes)
        byOutcome = ismember(SessionData.TrialOutcome(1:SessionData.nTrials), outcomes);
        out = out & byOutcome;
    end
    out = find(out);