%% Animal Specific Blocks
% DAT_9: April 18:22 25:28     (2 days where animal weak and not behaving
% well omitted)
%% DAT_(
saveName = 'DAT_9_TE';
sessionPath = 'E:\Data\DAT_9\SO_RewardPunish_odor\Session Data';
sessionNames = {...
    'DAT_9_SO_RewardPunish_odor_Apr16_2016_Session1',...
    'DAT_9_SO_RewardPunish_odor_Apr18_2016_Session2',...
    'DAT_9_SO_RewardPunish_odor_Apr19_2016_Session3',...
    'DAT_9_SO_RewardPunish_odor_Apr20_2016_Session1',...
    'DAT_9_SO_RewardPunish_odor_Apr21_2016_Session1',...
    'DAT_9_SO_RewardPunish_odor_Apr22_2016_Session1',...
    'DAT_9_SO_RewardPunish_odor_Apr25_2016_Session1',...
    'DAT_9_SO_RewardPunish_odor_Apr26_2016_Session1',...
    'DAT_9_SO_RewardPunish_odor_Apr27_2016_Session1',...
    'DAT_9_SO_RewardPunish_odor_Apr28_2016_Session1'...
    };

%% ChAT_16:
saveName = 'ChAT_16_TE';
sessionPath = 'E:\Data\ChAT_16\SO_RewardPunish_odor\Session Data';
sessionNames = {...
    'ChAT_16_SO_RewardPunish_odor_Apr16_2016_Session1.mat',...
    'ChAT_16_SO_RewardPunish_odor_Apr17_2016_Session1.mat',...
    'ChAT_16_SO_RewardPunish_odor_Apr18_2016_Session1.mat',...
    'ChAT_16_SO_RewardPunish_odor_Apr20_2016_Session1.mat',...
    'ChAT_16_SO_RewardPunish_odor_Apr21_2016_Session1.mat',...
    'ChAT_16_SO_RewardPunish_odor_Apr22_2016_Session1.mat',...
    'ChAT_16_SO_RewardPunish_odor_Apr23_2016_Session1.mat',...
    'ChAT_16_SO_RewardPunish_odor_Apr25_2016_Session1.mat',...    
    'ChAT_16_SO_RewardPunish_odor_Apr26_2016_Session1.mat',...    
    'ChAT_16_SO_RewardPunish_odor_Apr27_2016_Session1.mat',...    
    'ChAT_16_SO_RewardPunish_odor_Apr28_2016_Session1.mat'...        
    };


%% DAT_2
saveName = 'DAT_2_TE';
sessionPath = 'E:\Data\DAT_2\SO_RewardPunish_odor\Session Data';
sessionNames = {...
    'DAT_2_SO_RewardPunish_odor_Apr03_2016_Session2.mat',...
    'DAT_2_SO_RewardPunish_odor_Apr04_2016_Session2.mat',...
    'DAT_2_SO_RewardPunish_odor_Apr05_2016_Session1.mat',...
    'DAT_2_SO_RewardPunish_odor_Apr06_2016_Session1.mat'...    
    };


%%

clear TE;
TE = struct();
%%
for pos = 1:length(sessionNames);
    session = bpLoadSession(sessionNames{pos}, sessionPath);
    session.SessionData = demodulateSession(session.SessionData);
    nTrials = session.SessionData.nTrials;

    TE(pos).TrialTypes = session.SessionData.TrialTypes;
    TE(pos).TrialOutcome = session.SessionData.TrialOutcome;
    TE(pos).blLicks = bpCountEventByStates(session, 'Port1In', 'Cue', 'window', [-2 0]); % 2 seconds prior to cue
    TE(pos).anticipatoryLicks1 =  bpCountEventByStates(session, 'Port1In', 'Cue', 'endState', 'Delay'); % all licks cue + delay
    TE(pos).anticipatoryLicks2 = bpCountEventByStates(session, 'Port1In','Delay', 'referenceFromEnd', 1, 'window', [-2 0]);    % fixed 2 second window prceding US
    TE(pos).cueLicks = bpCountEventByStates(session, 'Port1In', 'Cue', 'window', [0.4 1.4]); % cue licks, matches phCuePeak
    TE(pos).cueLicksLong = bpCountEventByStates(session, 'Port1In', 'Cue', 'window', [0.4 2.4]); % cue licks, matches phCuePeak    
    TE(pos).delayLicks = bpCountEventByStates(session, 'Port1In', 'Delay', 'window', [0 1]); % delay licks (just 1st second), matches phDelayLicks, note that delay is somewhat variable
    TE(pos).usLicks = bpCountEventByStates(session, 'Port1In', 'Delay', 'referenceFromEnd', 1, 'window', [0 1]); % us licks, matches phUsPeak
    TE(pos).Photometry = processTrialAnalysis_Photometry(session);
    TE(pos).states = bpAddStatesAsTrialEvents(session);
    if isfield(session.SessionData, 'Epoch');
        TE(pos).epoch = session.SessionData.Epoch;
    else
        TE(pos).epoch = ones(1, nTrials);
    end
    TE(pos).fileName = repmat({session.filename}, nTrials, 1);
    TE(pos).phCuePeak = bpCalcPeak_dFF(TE(pos).Photometry, 1, [0.4 1.4], TE(pos).states.Cue); % cue is always 1s
    TE(pos).phCuePeakLong = bpCalcPeak_dFF(TE(pos).Photometry, 1, [0.4 2.4], TE(pos).states.Cue); % cue is always 1s
    TE(pos).phDelayPeak = bpCalcPeak_dFF(TE(pos).Photometry, 1, [0 1], TE(pos).states.Delay); % my delay is variable (1s initially, 2s for most), still I'm just taking 1st second
    TE(pos).phUsPeak = bpCalcPeak_dFF(TE(pos).Photometry, 1, [0 1], TE(pos).states.Delay, 'referenceFromEnd', 1);
    TE(pos).phBaseline = bpCalcPeak_dFF(TE(pos).Photometry, 1, [-2 0], TE(pos).states.Cue);    % baseline 2s before cue
end

%%
cd(sessionPath);
save(saveName, 'TE');

%%
allTE = struct(...
    'TrialTypes', [],...
    'TrialOutcome', [],...
    'blLicks', [],...
    'anticipatoryLicks1', [],...
    'anticipatoryLicks2', [],...
    'cueLicks', [],...
    'cueLicksLong', [],...
    'delayLicks', [],...
    'usLicks', [],...
    'epoch', [],...
    'fileName', {},...
    'phCuePeak', [],...
    'phCuePeakLong', [],...
    'phDelayPeak', [],...
    'phUsPeak', [],...
    'phBaseline', []...
    );
%     'dFF', {},... % ch1

for counter = 1:length(TE)
    allTE(1).TrialTypes = [allTE.TrialTypes ; TE(counter).TrialTypes(1:length(TE(counter).TrialOutcome))'];
    allTE(1).TrialOutcome = [allTE.TrialOutcome; TE(counter).TrialOutcome'];
    allTE(1).blLicks = [allTE.blLicks; TE(counter).blLicks.rate'];
    allTE(1).anticipatoryLicks1 = [allTE.anticipatoryLicks1; TE(counter).anticipatoryLicks1.rate'];
    allTE(1).anticipatoryLicks2 = [allTE.anticipatoryLicks2; TE(counter).anticipatoryLicks2.rate'];
    allTE(1).cueLicks = [allTE.cueLicks; TE(counter).cueLicks.rate'];
    allTE(1).cueLicksLong = [allTE.cueLicksLong; TE(counter).cueLicksLong.rate'];    
    allTE(1).delayLicks = [allTE.delayLicks; TE(counter).delayLicks.rate'];
    allTE(1).usLicks = [allTE.usLicks; TE(counter).usLicks.rate'];
    allTE(1).epoch = [allTE.epoch ; TE(counter).epoch'];
    allTE(1).fileName = [allTE.fileName ; TE(counter).fileName];
    allTE(1).phCuePeak = [allTE.phCuePeak ; TE(counter).phCuePeak.data];
    allTE(1).phCuePeakLong = [allTE.phCuePeakLong ; TE(counter).phCuePeakLong.data];    
    allTE(1).phDelayPeak = [allTE.phDelayPeak ; TE(counter).phDelayPeak.data];    
    allTE(1).phUsPeak = [allTE.phUsPeak ; TE(counter).phUsPeak.data];  
    allTE(1).phBaseline = [allTE.phBaseline ; TE(counter).phBaseline.data];  
end
%%
% convert fileName into indices
allNames= unique(allTE.fileName);
sessionIndex = zeros(size(allTE.TrialOutcome)); % now unique names are indices to seprates sessions
for counter = 1:length(allNames)
    sessionIndex(strcmp(allNames{counter}, allTE.fileName)) = counter;
end
allTE.sessionIndex = sessionIndex;
allTE.trialIndex = 1:length(allTE.epoch);

plotFields = {'blLicks', 'anticipatoryLicks1', 'anticipatoryLicks2', 'usLicks', 'epoch', 'sessionIndex',...
    'phCuePeak', 'phCuePeakLong', 'phDelayPeak', 'phUsPeak', 'trialIndex'};
clear cuedReward 
clear cuedPunish
cuedReward.Trials = bpFilterTrials2(allTE, 'TrialTypes', 1);
cuedPunish.Trials = bpFilterTrials2(allTE, 'TrialTypes', 3);
for counter = 1:length(plotFields)
    cuedReward.(plotFields{counter}) = allTE.(plotFields{counter})(cuedReward.Trials);
    cuedPunish.(plotFields{counter}) = allTE.(plotFields{counter})(cuedPunish.Trials);    
end


%%
ensureFigure([saveName '_longitudinalCueResponses'], 1); 
subplot(4,1,1); scatter(cuedReward.trialIndex, cuedReward.anticipatoryLicks2 - cuedReward.blLicks,'.'); % both blLicks and anticipatoryLicks2 are from 2s periods
% subplot(4,1,1); scatter(1:length(allTE.epoch), cuedReward.anticipatoryLicks2);
maxLR = max(cuedReward.anticipatoryLicks2);
% subplot(4,1,1); plot(smooth(cuedReward.anticipatoryLicks2 / mean(cuedReward.blLicks)));
% maxLR = smooth(max(cuedReward.anticipatoryLicks2 / mean(cuedReward.blLicks)));
hold on; stem(cuedReward.trialIndex(1:end-1), abs(diff(cuedReward.epoch)) * maxLR, 'g', 'Marker', 'none');
hold on; stem(cuedReward.trialIndex(1:end-1), diff(cuedReward.sessionIndex) * maxLR, 'r', 'Marker', 'none');
ylabel('antic. licks'); title('cued reward trials');
set(gca, 'XLim', [1 length(allTE.epoch)]);

subplot(4,1,2); scatter(cuedReward.trialIndex, cuedReward.phCuePeak, '.');
maxP = max(cuedReward.phCuePeak);
hold on; stem(cuedReward.trialIndex(1:end-1), abs(diff(cuedReward.epoch)) * maxP, 'g', 'Marker', 'none');
hold on; stem(cuedReward.trialIndex(1:end-1), diff(cuedReward.sessionIndex) * maxP, 'r', 'Marker', 'none');
ylabel('cue dF/F');
set(gca, 'XLim', [1 length(allTE.epoch)]); set(gca, 'YLim', [-0.2 0.2]);



subplot(4,1,3); scatter(cuedPunish.trialIndex, cuedPunish.anticipatoryLicks2 - cuedPunish.blLicks, '.');
maxLR = max(cuedPunish.anticipatoryLicks2);
hold on; stem(cuedPunish.trialIndex(1:end-1), abs(diff(cuedPunish.epoch)) * maxLR, 'g', 'Marker', 'none');
hold on; stem(cuedPunish.trialIndex(1:end-1), diff(cuedPunish.sessionIndex) * maxLR, 'r', 'Marker', 'none');
ylabel('antic. licks'); title('cued punish trials');
set(gca, 'XLim', [1 length(allTE.epoch)]);

subplot(4,1,4); scatter(cuedPunish.trialIndex, cuedPunish.phCuePeak, '.');
maxP = min(cuedPunish.phCuePeak);
hold on; stem(cuedPunish.trialIndex(1:end-1), abs(diff(cuedPunish.epoch)) * maxP, 'g', 'Marker', 'none');
hold on; stem(cuedPunish.trialIndex(1:end-1), diff(cuedPunish.sessionIndex) * maxP, 'r', 'Marker', 'none');
ylabel('cue dF/F');
saveas(gcf, [saveName '_longitudinalCueResponses.fig']);
disp('figure saved');
set(gca, 'XLim', [1 length(allTE.epoch)]);


%%
cuedreward10 = bpFilterTrials2(TE(10), 'TrialTypes', 1);
ensureFigure('DAT_9_April_28', 1); scatter(TE(10).phCuePeak.data(cuedreward10), TE(10).phUsPeak.data(cuedreward10));

%% for DAT_2 trial by trial

plotFields = {'blLicks', 'anticipatoryLicks1', 'anticipatoryLicks2', 'usLicks', 'epoch', 'sessionIndex',...
    'phCuePeak', 'phCuePeakLong', 'phDelayPeak', 'phUsPeak', 'trialIndex'};
clear cuedReward 
cuedReward.Trials = bpFilterTrials2(allTE, 'TrialTypes', 1, 'fileName', 'DAT_2_SO_RewardPunish_odor_Apr03_2016_Session2.mat');
for counter = 1:length(plotFields)
    cuedReward.(plotFields{counter}) = allTE.(plotFields{counter})(cuedReward.Trials);
end

ensureFigure('DAT_2_trialByTrial', 1);
subplot(2,2,1);
scatter(cuedReward.anticipatoryLicks2, cuedReward.phCuePeak); xlabel('Antic. Licks'); ylabel('Cue, dF/F');
subplot(2,2,2);
scatter(cuedReward.anticipatoryLicks2, cuedReward.phUsPeak); xlabel('Antic. Licks'); ylabel('Reward, dF/F');
subplot(2,2,3);
scatter(cuedReward.phCuePeak, cuedReward.phUsPeak); xlabel('Cue, dF/F'); ylabel('Reward, dF/F');
subplot(2,2,4);
scatter(cuedReward.phCuePeak, cuedReward.phUsPeak); xlabel('Cue, dF/F'); ylabel('Reward, dF/F'); hold on;
fitobject = fit(cuedReward.phCuePeak, cuedReward.phUsPeak, 'poly1');
xlabel('Cue, dF/F'); ylabel('Reward, dF/F'); 
plot(fitobject); title('Same as Left with linear fit');
saveas(gcf, 'DAT_2_trialByTrial.fig');
saveas(gcf, 'DAT_2_trialByTrial.tif');

%%
%% for DAT_2 trial by trial (included uncued trials to show that I get a negative correlation now)

plotFields = {'blLicks', 'anticipatoryLicks1', 'anticipatoryLicks2', 'usLicks', 'epoch', 'sessionIndex',...
    'phCuePeak', 'phCuePeakLong', 'phDelayPeak', 'phUsPeak', 'trialIndex'};
clear cuedReward 
cuedReward.Trials = bpFilterTrials2(allTE, 'TrialTypes', [1 5], 'fileName', 'DAT_2_SO_RewardPunish_odor_Apr03_2016_Session2.mat');
for counter = 1:length(plotFields)
    cuedReward.(plotFields{counter}) = allTE.(plotFields{counter})(cuedReward.Trials);
end

ensureFigure('DAT_2_trialByTrial', 1);
subplot(2,2,1);
scatter(cuedReward.anticipatoryLicks2, cuedReward.phCuePeak); xlabel('Antic. Licks'); ylabel('Cue, dF/F');
title('Cued and UNCUED reward trials');
subplot(2,2,2);
scatter(cuedReward.anticipatoryLicks2, cuedReward.phUsPeak); xlabel('Antic. Licks'); ylabel('Reward, dF/F');
subplot(2,2,3);
scatter(cuedReward.phCuePeak, cuedReward.phUsPeak); xlabel('Cue, dF/F'); ylabel('Reward, dF/F'); hold on;
subplot(2,2,4);
scatter(cuedReward.phCuePeak, cuedReward.phUsPeak); xlabel('Cue, dF/F'); ylabel('Reward, dF/F'); hold on;
fitobject = fit(cuedReward.phCuePeak, cuedReward.phUsPeak, 'poly1');
xlabel('Cue, dF/F'); ylabel('Reward, dF/F'); 
plot(fitobject); title('Same as Left with linear fit');
saveas(gcf, 'DAT_2_trialByTrial_plusUncued.fig');
saveas(gcf, 'DAT_2_trialByTrial_plusUncued.tif');