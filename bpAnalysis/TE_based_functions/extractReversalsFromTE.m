function split = extractReversalsFromTE(TE, trials, inputData, varargin)
% TE- TrialEvents structure that contains fields BlockChange and filename
% trials: logical or linear indices
% varargin: name (string), value (nGlobalTrials x 1) pairs corresponding to
% fields in TE that you want to parse and deposit into a reversal split
% output.  example: 
% TE.reversals = extractReversalsFromTE(TE, rewardTrials & hitTrials,...
% 'csLicks', TE.csLicks.rate, 'csDFF_ch1', TE.phPeakMean_cs(1).data);
% TE.reversals has fields 'csLicks', 'csDFF_ch1', and 'globalTrialNumber'
% Each field has subfields .before and .after with 1 row for each reversal
% (.before contains data from pre-reversal block of trials, and .after from
% post-reversal block of trials)

% (probably a more efficient way to do this, annoying mixed use of logical
% and linear indexing)

    defaults = {...
        'maxReversals', 1;... % maximum number of reversals to consider per session
        };

    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings

    split = struct();    
    
    emptyTrials = false(size(TE.filename));
    if ~islogical(trials)
        trials2 = emptyTrials;
        trials2(trials) = 1;
        trials = true(trials2);
    end
    
    % find reversal and session breaks
    % get start and stop indices of before and after trials for each
    % reversal, and largest pre- and post- reversal trial blocks
    maxTrials = zeros(1,2);
    rev = find(TE.BlockChange);
    before = zeros(1, 2);
    after = zeros(1, 2);
    si = [1; find(TE.sessionChange); length(TE.filename) + 1];
    revCounter = 1;
    
    for scounter = 1:length(si) - 1
        s1 = si(scounter);
        s2 = si(scounter + 1);
        % 4 indices associated with each reversal STARTbefore, ENDbefore,
        % STARTafter, ENDafter
        theseReversals = [s1; rev(rev > s1 & rev < s2); s2]; % NOT >= because you don't want reversals occuring at beginning or end of session (not that these should occur)
        for counter = 2:length(theseReversals) - 1
            before(revCounter,1) = theseReversals(counter - 1);
            before(revCounter, 2) = theseReversals(counter) - 1;
            after(revCounter,1) = theseReversals(counter);
            after(revCounter,2) = theseReversals(counter + 1) -1;
            maxTrials(1) = max(maxTrials(1), sum(trials(before(1):before(2)))); % max before block            
            maxTrials(2) = max(maxTrials(1), sum(trials(after(1):after(2)))); % max after block                        
            revCounter = revCounter + 1;
            if counter == s.maxReversals + 1
                break
            end
        end
    end
    nRev = size(before, 1);
    % theseTrials = zeros(size(trials));
    % theseTrials(before(1):before(2) = trials;

    %% compile reversal blocks for desired trial events    
    inputData = [{'globalTrialNumber', 1:length(TE.filename)} inputData];
    counter = 1;
    
    while counter+1 <= length(inputData) 
        fieldName = inputData{counter};
        te = inputData{counter+1}; 
        % initialize
        if iscell(te)
            split.(fieldName).before = cell(nRev, maxTrials(1));
            split.(fieldName).after = cell(nRev, maxTrials(2));            
        else
            split.(fieldName).before = NaN(nRev, maxTrials(1));
            split.(fieldName).after = NaN(nRev, maxTrials(2));
        end
        % compile before and after blocks
        for revCounter = 1:nRev
            % before
            theseTrials = emptyTrials;
            theseTrials(before(revCounter, 1):before(revCounter,2)) = trials(before(revCounter, 1):before(revCounter,2));
            revData = te(theseTrials);
            split.(fieldName).before(revCounter, end - length(revData) + 1:end) = revData;
            % after
            theseTrials = emptyTrials;
            theseTrials(after(revCounter, 1):after(revCounter,2)) = trials(after(revCounter, 1):after(revCounter,2));
            revData = te(theseTrials);
            split.(fieldName).after(revCounter, 1:length(revData)) = revData;      
        end
        counter=counter+2;
    end
    split.trialsBefore = linspace(-maxTrials(1)+1, 0, maxTrials(1));
    split.trialsAfter = linspace(1,maxTrials(2),maxTrials(2));
    