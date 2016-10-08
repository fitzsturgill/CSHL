 
function TE = truncateSessionsFromTE(TE, action)
    global TRUNC
    if nargin < 2
        tei = 1;
    end
    
    switch action
        case 'init'
            nSessions = max(TE.sessionIndex);
            nTrials = length(TE.sessionIndex);            
            lastTrial = [find(TE.sessionChange) - 1, ; nTrials];
            firstTrial = [0 find(TE.sessionChange)];
            TRUNC = structure(...
                'nSessions', nSessions,...
                'lastTrial', lastTrial,...
                'firstTrial', firstTrial,...
                'truncTrial', lastTrial,... % begin with trunc trial indicators pointing at last trial in each session
                'fig', [],...
                'ax', [],...
                'licksHandle', [],... % handle for reward licks vs trial number line plot
                'lastTrialHandle', [],... % handle for trunc trial indicator
                'currentSession', 1, ... % index of session being interactively adjusted
                'reject', zeros(1, nTrials)...
                );
                

            %% plot reward period lick rate vs trial number to visualize satiation/ lapsing behavior towards end of each session
            TRUNC.fig = ensureFigure('truncFigure', 1);
            TRUNC.rewardTrials = find(filterTE(TE, 'trialOutcome', 1));
            TRUNC.rewardLicks = smooth(TE(tei).usLicks.rate(rewardTrials), 5);
            
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
                end
            end
            plot(find(rewardTrials), smooth(TE(tei).usLicks.rate(rewardTrials), 5)); hold on; 
            plot(find(rewardTrials), [0; diff(TE(tei).sessionIndex(rewardTrials))] * 10);
            ylabel('Lick/s, Us'); xlabel('trial #'); textBox(TE(tei).filename{1}(1:7));
            set(gca, 'YLim', [0 10]);                
        case 'update'
            TE.reject = TRUNC.reject;
    end
end







