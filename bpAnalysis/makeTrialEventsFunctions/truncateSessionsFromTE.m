 
function TE = truncateSessionsFromTE(TE, action)
    % interactively adjust session truncation points
    % designed with cuedOutcome_Odor_Complete TE structure in mind
    % needs to be updated for other TEs- i.e. don't hardcode rewardLicks
    % TE field (currently usLicks)
    % WARNING! currently there are two ways to update
    % 1) call with 'update' as action, this will not overwrite TE, except
    % via output argument
    % 2) 'u' keypres- this WILL overwrite TE
    global TRUNC
    evalin('base', 'global TRUNC');
    % trunc trial is LAST INCLUDED TRIAL in each session
    yMax = 10;
    switch action
        case 'init'
            nSessions = max(TE.sessionIndex);
            nTrials = length(TE.sessionIndex);            
            lastTrial = [find(TE.sessionChange) - 1, ; nTrials]; % the last trial in each session
            firstTrial = [1; find(TE.sessionChange)]; % the first trial in each session
            TRUNC = struct(...
                'nSessions', nSessions,...
                'lastTrial', lastTrial,...
                'firstTrial', firstTrial,...
                'truncTrial', lastTrial,... % begin with trunc trial indicators pointing at last trial in each session
                'fig', [],...
                'ax', [],...
                'licksHandle', [],... % handle for reward licks vs trial number line plot
                'truncTrialHandle', [],... 
                'currentSession', 1, ... % index of session being interactively adjusted
                'reject', zeros(1, nTrials)...
                );
                

            %% plot reward period lick rate vs trial number to visualize satiation/ lapsing behavior towards end of each session
            TRUNC.fig = ensureFigure('truncFigure', 1);
            set(TRUNC.fig, 'WindowKeyPressFcn', @truncKeyPressFcn,...
                'CloseReq', @truncClose);
            TRUNC.rewardTrials = find(filterTE(TE, 'trialOutcome', 1));
            TRUNC.rewardLicks = smooth(TE.usLicks.rate(TRUNC.rewardTrials), 5);
            TRUNC.truncTrialHandle = zeros(1, nSessions); % will contain handles for trunc trial indicators
            
            TRUNC.licksHandle = plot(TRUNC.rewardTrials, smooth(TE.usLicks.rate(TRUNC.rewardTrials), 5)); hold on; 
            plot(TRUNC.rewardTrials, [0; diff(TE.sessionIndex(TRUNC.rewardTrials))] * yMax);
            % plot trunc trial indicators
            for session = 1:TRUNC.nSessions
                truncIndex = nearest(TRUNC.rewardTrials, TRUNC.truncTrial(session)); % take nearest reward outcome trial to trunc trial, better to step along along trials than just reward trials
                TRUNC.truncTrialHandle(session) = line('XData', TRUNC.rewardTrials(truncIndex), 'YData', TRUNC.rewardLicks(truncIndex), 'Marker', 'o',...
                    'MarkerSize', 8,...
                    'MarkerFacecolor', 'm'); 
            end
            set(TRUNC.truncTrialHandle(TRUNC.currentSession), 'MarkerFaceColor', 'g'); % highlight current session trunc marker
                
            ylabel('Lick/s, Us'); xlabel('trial #'); textBox(TE.filename{1}(1:7));
            set(gca, 'YLim', [0 yMax]);
            
            % attempt to use truncation points implied by reject field in
            % TE (if it exists)
            if isfield(TE, 'reject')
                TRUNC.reject = TE.reject;
                for session = 1:TRUNC.nSessions
                    counter = TRUNC.lastTrial(session);
                    % scan back from last trial and find last rejected
                    % trial (might be able to do this without a loop but whatever)
                    while TRUNC.reject(counter) == 1
                        if counter == TRUNC.firstTrial(session);
                            break % if for some reason all trials in session have been rejected, stop at first trial
                        end
                        counter = counter - 1;
                    end
                    TRUNC.truncTrial(session) = counter;
                    updateTrunc(session, counter); % update
                end
            else
                TE.reject = TRUNC.reject;
            end

        case 'update'
            TE.reject = TRUNC.reject;
    end
end

function truncKeyPressFcn(src, evt) % window key press fcn, executes whenever figure or its children has/have focus...
    global TRUNC
    si = TRUNC.currentSession; % session index
    switch evt.Character
        case 28 % left arrow
            if ~ismember('shift', evt.Modifier)
                updateTrunc(si, TRUNC.truncTrial(si) - 1);
            else
                updateTrunc(si, TRUNC.truncTrial(si) - 10);
            end
        case 29 % right arrow
            if ~ismember('shift', evt.Modifier)
                updateTrunc(si, TRUNC.truncTrial(si) + 1);
            else
                updateTrunc(si, TRUNC.truncTrial(si) + 10);
            end   
        case 30 % up arrow
            if si + 1 <= TRUNC.nSessions
                updateTrunc(si + 1);
            else
                updateTrunc(1); % cycle to first session
            end
        case 31 % down arrow
            if si - 1 >= 1
                updateTrunc(si - 1);
            else
                updateTrunc(TRUNC.nSessions); % cycle to last session
            end
        case 117 % u, update 
            evalin('base', 'TE.reject=TRUNC.reject;');
            display('*** truncateSessionsFromTE: updated TE.reject ***');
        otherwise
    end
end

function updateTrunc(s,t)
    % s --> session number, pass [] to use current session
    % t --> new trial number of trunc trial for current session
    global TRUNC
    if nargin < 2
        t = 0;
    end
    if isempty(s)
        s = TRUNC.currentSession;
    end
    if TRUNC.currentSession ~= s
        set(TRUNC.truncTrialHandle(TRUNC.currentSession), 'MarkerFaceColor', 'm'); % unhighlight old current session trunc marker
        TRUNC.currentSession = s;
        set(TRUNC.truncTrialHandle(TRUNC.currentSession), 'MarkerFaceColor', 'g'); % highlight current session trunc marker
    end
    
    if t
        TRUNC.truncTrial(s) = min(max(t, TRUNC.firstTrial(s)), TRUNC.lastTrial(s));
        truncIndex = nearest(TRUNC.rewardTrials, TRUNC.truncTrial(s));
        set(TRUNC.truncTrialHandle(s), 'XData', TRUNC.rewardTrials(truncIndex), 'YData', TRUNC.rewardLicks(truncIndex));
        
        % update reject field, trunc trial is LAST INCLUDED TRIAL        
        TRUNC.reject(TRUNC.firstTrial(s):TRUNC.truncTrial(s)) = 0;
        if TRUNC.truncTrial(s) < TRUNC.lastTrial(s)
            TRUNC.reject((TRUNC.truncTrial(s) + 1):TRUNC.lastTrial(s)) = 1;
        end
    end
end
    
    
    
function truncClose(src,evt)
    global TRUNC
    clf(TRUNC.fig);
    delete(TRUNC.fig);
    clear TRUNC;
end




