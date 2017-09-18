function TE = makeTE_CuedOutcome_Odor_Complete_Nlx(sessionPath, sessionName)
    
    if nargin < 2
        if nargin == 1
            cd(sessionPath);
        end
        sessions = bpLoadSessions; % select Bpod session file from within CellBase file structure, e.g. within C:\FitzData\Cellbase_dev\CD4\170909a
    else
        sessions = bpLoadSessions([], sessionName, sessionPath);
    end
    if isempty(sessions)
        return
    end
    if length(sessions) > 1
        error('*** only handles single sessions for now ***');
    end
    TE = makeTE_CuedOutcome_Odor_Complete(sessions, 1);
    
    load(fullfile(sessions.filepath, 'Events.mat'));   % load Neuralynx events
    % now Events_EventIDs, Events_EventStrings, Events_Extras,
    % Events_Nttls, and Events_TimeStamps are loaded into the local
    % scope/workspace
    
    if nargin < 1
        sessionPath = sessions.filepath;
    end

    
    
    TrialStart_nlx = getBehaviorStartTimes(Events_Nttls, Events_EventStrings, Events_TimeStamps);
    TrialStart_Bpod = sessions.SessionData.TrialStartTimestamp;
    % Match timestamps - in case of mismatch, try to fix
if ~ismatch(TrialStart_Bpod,TrialStart_nlx)
    % note: obsolete due the introduction of TTL parsing
    TrialStart_nlx = clearttls(TrialStart_nlx); % eliminate recorded TTL's within 0.5s from each other - broken TTL pulse
    if ~ismatch(TrialStart_Bpod,TrialStart_nlx)
        TrialStart_nlx = trytomatch(TrialStart_Bpod,TrialStart_nlx);  % try to match time series by shifting
        if ~ismatch(TrialStart_Bpod,TrialStart_nlx)
            TrialStart_nlx = tryinterp(TrialStart_Bpod,TrialStart_nlx); % interpolate missing TTL's or delete superfluous TTL's up to 10 erroneous TTl's
            if ~ismatch(TrialStart_Bpod,TrialStart_nlx)  % TTL matching failure
                error('MakeTrialEvents:TTLmatch','Matching TTLs failed.')
            else
                warning('MakeTrialEvents:TTLmatch','Missing TTL interpolated.')
            end
        else
            warning('MakeTrialEvents:TTLmatch','Shifted TTL series.')
        end
    else
        warning('MakeTrialEvents:TTLmatch','Broken TTLs cleared.')
    end
end

% Eliminate last TTLs recorded in only one system (Bpod or Nlx)
if length(TrialStart_nlx) > length(TrialStart_Bpod)
    TrialStart_nlx = TrialStart_nlx(1:length(TrialStart_Bpod));
elseif length(TrialStart_nlx) < length(TrialStart_Bpod)
    TE = shortenTE(TE, length(TrialStart_nlx));
end

 TE.TrialStart = TrialStart_nlx(:);
     
%     if sessions.SessionData.nTrials ~= length(TE.TrialStart)
%         error('*** mismatched numbers of trials: in makeTE_CuedOutcome_Odor_Complete_Nlx ***');
%     end
    
    save(fullfile(sessionPath, 'TrialEvents.mat'), 'TE')
    disp(['*** Saved: ' fullfile(sessionPath, 'TrialEvents.mat') ' ****']);
    
    
    
    
    
    function I = ismatch(ts,TrialStart_nlx)

% Check if the two time series match notwithstanding a constant drift
clen = min(length(ts),length(TrialStart_nlx));
I = abs(max(diff(ts(1:clen)-TrialStart_nlx(1:clen)))) < 0.1;  % the difference between the timestamps on 2 systems may have a constant drift, but it's derivative should still be ~0

% note: abs o max is OK, the derivative is usually a small neg. number due
% to drift of the timestamps; max o abs would require a higher tolerance
% taking the drift into account (if 2 event time stamps are far, the drift
% between them can be large)

% -------------------------------------------------------------------------
function TrialStart_nlx = tryinterp(ts,TrialStart_nlx)

% Interpolate missing TTL's or delete superfluous TTL's up to 10 erroneous TTl's
for k = 1:10
    if ~ismatch(ts,TrialStart_nlx)
        son3 = TrialStart_nlx - TrialStart_nlx(1) + ts(1);
        adt = diff(ts(1:min(length(ts),length(TrialStart_nlx)))-TrialStart_nlx(1:min(length(ts),length(TrialStart_nlx))));
        badinx = find(abs(adt)>0.1,1,'first') + 1;  % find problematic index
        if adt(badinx-1) < 0    % interploate
            ins = ts(badinx) - linterp([ts(badinx-1) ts(badinx+1)],[ts(badinx-1)-son3(badinx-1) ts(badinx+1)-son3(badinx)],ts(badinx));
            TrialStart_nlx = [TrialStart_nlx(1:badinx-1) ins+TrialStart_nlx(1)-ts(1) TrialStart_nlx(badinx:end)];
        else
%             ins = son3(badinx) - linterp([son3(badinx-1) son3(badinx+1)],[son3(badinx-1)-ts(badinx-1) son3(badinx+1)-ts(badinx)],son3(badinx));
%             ts = [ts(1:badinx-1) ins ts(badinx:end)];
            TrialStart_nlx(badinx) = [];   % delete
        end
    end
end

% -------------------------------------------------------------------------
function TrialStart_nlx = trytomatch(ts,TrialStart_nlx)

% Try to match time series by shifting
len = length(TrialStart_nlx) - 15;
minx = nan(1,len);
for k = 1:len
    minx(k) = max(diff(ts(1:15)-TrialStart_nlx(k:k+14)));  % calculate difference in the function of shift
end
mn = min(abs(minx));
minx2 = find(abs(minx)==mn);
minx2 = minx2(1);   % find minimal difference = optimal shift
TrialStart_nlx = TrialStart_nlx(minx2:min(minx2+length(ts)-1,length(TrialStart_nlx)));

% -------------------------------------------------------------------------
function TrialStart_nlx = clearttls(TrialStart_nlx)

% Eliminate recorded TTL's within 0.5s from each other
inx = [];
for k = 1:length(TrialStart_nlx)-1
    s1 = TrialStart_nlx(k);
    s2 = TrialStart_nlx(k+1);
    if s2 - s1 < 0.5
        inx = [inx k+1]; %#ok<AGROW>
    end
end
TrialStart_nlx(inx) = [];

% -------------------------------------------------------------------------
function TE2 = shortenTE(TE2,n)

% Eliminate behavioral trials
fnm = fieldnames(TE2);
for k = 1:length(fieldnames(TE2))
    TE2.(fnm{k}) = TE2.(fnm{k})(1:n);
end