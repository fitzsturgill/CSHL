function [lickTimes, lickTrials, Cavg, Cx, Csem] = researchStatementLicks(sessions, type, outcome)
% crapppy, for research statemntn

    
    lickTimes = [];
    lickTrials = [];
    
    
    % histogram from 1sec to 5seconds
    deltaT = 0.5;
    Cx = -2:deltaT:7;
    C = [];
    for j = 1:length(sessions)
        session = sessions(j).SessionData;
        trials = bpFilterTrials(session, type, outcome);
        length(trials)
        for i = 1:length(trials)
            trial = trials(i);
            if isfield(session.RawEvents.Trial{trial}.Events, 'Port2In')
                theseLicks = session.RawEvents.Trial{trial}.Events.Port2In;
                stimStart = session.RawEvents.Trial{trial}.States.DeliverStimulus(1);
                theseLicks = theseLicks - stimStart;
                lickTimes = [lickTimes theseLicks];
                if j == 1 && i == 1
                    lickTrials = [lickTrials zeros(1, length(theseLicks)) + 1];
                else
                    lickTrials = [lickTrials zeros(1, length(theseLicks)) + lickTrials(end) + 1];
                end

            else % there are no licks
                theseLicks = NaN;
                if j == 1 && i == 1
                    lickTimes(1) = theseLicks; % 
                    lickTrials(1) = 1;
                else
                    lickTimes(end + 1) = theseLicks;
                    lickTrials(end + 1) = lickTrials(end) + 1;
                end
                
            end
            if j == 1 && i == 1
                C = histcounts(theseLicks, Cx);
%                 C = smooth(C')';
%                 try
%                     C = smooth(C, 3);
%                 end
            else
                newline = histcounts(theseLicks, Cx);
%                 newline = smooth(newline')';
%                 try
%                     newl = smooth(newl, 3);
%                 catch
%                     disp('wtf');
%                 end
                C = [C; newline];
            end
        
        end
        
    end
    %return the bin centers for plotting
    Cx = Cx(2:end) - diff(Cx(1:2));
    
    
    disp(['num of trials is :' num2str(size(C, 1))]);
    Crates = C;
    Crates = C ./ size(C, 1) ./ deltaT;
    Cavg = median(Crates);
    Csem = std(Crates) ./ sqrt(size(Crates, 1));
    
    
    
    
    
    
    