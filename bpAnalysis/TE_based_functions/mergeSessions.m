function TE = mergeSessions(TE, sessionIndicesToMerge)
    % sessionsToMerge-   consecutive sessions (as loaded in TE), denoted by sessionIndex  to merge

    sessionIndicesToMerge = sort(sessionIndicesToMerge);    
    newIndex = TE.sessionIndex(find(TE.sessionIndex == sessionIndicesToMerge(1), 1));
    newName = TE.filename(find(TE.sessionIndex == sessionIndicesToMerge(1), 1));
    TE.sessionIndex(ismember(TE.sessionIndex, sessionIndicesToMerge)) = newIndex;
    TE.filename(ismember(TE.sessionIndex, sessionIndicesToMerge)) = newName;
    if max(unique(TE.sessionIndex) > max(sessionIndicesToMerge)) % meaning you're not just merging the last sessions        
        shiftTheseDown = max(sessionIndicesToMerge) + 1 : max(TE.sessionIndex);
        changeToThis = (1:length(shiftTheseDown)) + newIndex;
        
        for counter = 1:length(shiftTheseDown)
            TE.sessionIndex(TE.sessionIndex == shiftTheseDown(counter)) = changeToThis(counter);
        end
    end
    
    TE.sessionChange = [0; diff(TE.sessionIndex)];
    