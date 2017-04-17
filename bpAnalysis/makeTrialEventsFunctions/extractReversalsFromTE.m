function split = extractReversalsFromTE(TE, trials, varargin)
    % trials: logical or linear indices

    
    split = structure();    
    
    % find reversal and session breaks
    if islogical(trials)
        trials = find(trials);
    end
    
    % get start and stop indices of before and after trials for each
    % reversal
    maxTrials = zeros(1,2);
    rev = find(TE.BlockChange);
    before = zeros(length(blockChanges), 2);
    after = zeros(length(blockChanges), 2);
    si = [1 find(TE.sessionChange) length(TE.filename) + 1];
    revCounter = 1;
    s1=0;s2=0;
    
    % EASY WAY TO FIND MAX BLOCKS- JUST ASSUME THAT MAX SIZE IS EQUAL TO
    % THE NUMBER OF TRIALS IN THE BIGGEST SESSION- IT'LL BE SMALLER BUT
    % WE'RE NOT DEALING WITH HUGE DATA ARRAYS (NOT MEMORY LIMITED). YOU CAN
    % UNPAD THE ARRAY AT THE END
    % 
    % HOW TO USE THE INDICES.   theseTrials = zeros(size(trials));
    % theseTrials(before(1):before(2) = trials;
    for scounter = 1:length(si) - 1
        s1 = si(scounter);
        s2 = si(scounter + 1);
        % 4 indices associated with each reversal STARTbefore, ENDbefore,
        % STARTafter, ENDafter
        theseReversals = [s1 rev(rev > s1 & rev < s2) s2]; % NOT >= because you don't want reversals occuring at beginning or end of session (not that these should occur)
        for counter = 2:length(theseReversals) - 1
            before(revCounter,1) = theseReversals(counter - 1);
            before(revCounter, 2) = theseReversals(counter) - 1;
            after(revCounter,1) = theseReversals(counter);
            after(revCounter,2) = theseReversals(counter + 1) -1;
            revCounter = revCounter + 1;
            
            maxTrials(1) = max(maxTrials(1), find(trials >= before(1) & trials <= before(2)));
        end
    end
    
    % find the biggest before and after reversal blocks
    
         
    
    %% STopping point- Now I have before, and after, 2 lists of start and stop points corresponding to before reversal epochs and after reversal epochs.
            
            
    
    
    
    
    % parse new field/ trial event parameter pairs    
    counter = 1;
    while counter+1 <= length(varargin) 
        fieldName = varargin{counter};
        te = varargin{counter+1}; % must be numerical to begin with
        split.(
        counter=counter+2;
    end
    