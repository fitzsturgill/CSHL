
fig = ensureFigure('FSChAT10_SO_Training_FS_Jan07_2016_Session1', 1);
skipAnalysis = 0;
conditions = {'toneA_reward', 'toneA_omit', 'toneA_punish', 'toneB_reward', 'toneB_omit', 'toneB_punish'};
trialTypes = [ 1 1 1 2 2 2];
outcomes = [2 3 4 2 3 4];
linespecs = {'g', 'k', 'r', 'g', 'k', 'r'};
if ~skipAnalysis
    
    demodData = struct(...
        'condition', '',...
        'trialType', [],...
        'trialOutcome', [],...
        'avg', [],...
        'avgX', [],...
        'avgSEM', []...
        );
    demodData = repmat(demodData, 1, length(conditions));
end

subplot(2,1,1); hold on    
ylabel('GCAMP6s'); xlabel('time(s)'); title('tone A');
for counter = 1:length(conditions)

    if ~skipAnalysis
        demodData(counter).condition = conditions{counter};
        demodData(counter).trialType = trialTypes(counter);
        demodData(counter).trialOutcome = outcomes(counter);    
        [avg, avgX, avgSEM] = demodSession_Old(session.SessionData, trialTypes(counter), outcomes(counter));
        demodData(counter).avg = avg;
        demodData(counter).avgX = avgX;
        demodData(counter).avgSEM = avgSEM;
    end
    disp('*** condition done ***');
    if counter == 4
        subplot(2,1,2); hold on % switch to different axes for toneB
        ylabel('GCAMP6s'); xlabel('time(s)'); title('tone B');        
    end
    boundedline(demodData(counter).avgX, demodData(counter).avg, demodData(counter).avgSEM, linespecs{counter}, gca, 'alpha');     
end





