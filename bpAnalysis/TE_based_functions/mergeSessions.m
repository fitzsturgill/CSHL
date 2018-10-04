function TE = mergeSessions(TE, sessionIndicesToMerge)
    % sessionsToMerge-   consecutive sessions (as loaded in TE), denoted by sessionIndex  to merge
    
    sessionIndicesToMerge = sessionIndicesToMerge(:)'; % ensure row vector
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
    
    % if sessions field exists, update indices to reflect change
    if isfield(TE, 'sessions')
        for index = sessionIndicesToMerge(2:end)
            TE.sessions(index).index = newIndex;
        end
        if index < length(TE.sessions)
            for trailing = index+1:length(TE.sessions)
                TE.sessions(trailing).index = TE.sessions(trailing - 1).index + 1;
            end
        end
        
    end
    