saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_LNL_odor_V2(sessions);

%% extract peak trial dFF responses to cues and reinforcement and lick counts
% csWindow = cellfun(@(x,y,z) [x(1) max(y(end), z(end))], TE.Cue, TE.AnswerLick, TE.AnswerNoLick); % max- to select either AnswerLick or AnswerNoLick timestamp (unused state contains NaN)
nTrials = length(TE.filename);
csWindow = zeros(nTrials, 2);
csWindow(:,2) = cellfun(@(x,y,z) max(x(end), y(end)) - z(1), TE.AnswerLick, TE.AnswerNoLick, TE.Cue); 
% max 1) to select AnswerNoLick time stampfor no lick trials (unused state contains NaN)
% 2) to select AnswerLick time stamp for lick trials (AnswerLick follows
% AnswerNoLick state)


TE.csLicks = countEventFromTE(TE, 'Port1In', csWindow, TE.Cue);

usWindow = [0 0.5];
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'


TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);

%% 
basepath = 'Z:\SummaryAnalyses\LNL_odor_v2_new\noPunish_noPhotometry\';
sep = strfind(TE.filename{1}, '_');
subjectName = TE.filename{1}(1:sep(2)-1);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);


%%
truncateSessionsFromTE(TE, 'init', 'usLicks', filterTE(TE, 'trialType', 1));
%%  generate trial lookups for different combinations of conditions
LNL_conditions;

%% calculate discriminability of csPlus from csMinus trials according to anticipatory licking
TE.auROC = calc_auROC_reversalsFromTE(TE, TE.csLicks.rate, csPlusTrials, csMinusTrials, 20, 1);

%% plot licks vs reversals scatter plot
smoothFactor = 5;
saveName = [subjectName '_longitudinal_Licking'];
h=ensureFigure(saveName, 1); 
mcLandscapeFigSetup(h);
subplot(2,1,1); hold on
thisscatter = scatter(trialCount(csPlusTrials), TE.csLicks.rate(csPlusTrials), 10, 'o'); % both blLicks and anticipatoryLicks2 are from 2s periods
plot(trialCount(csPlusTrials), smooth(TE.csLicks.rate(csPlusTrials), smoothFactor), 'b', 'LineWidth', 2);
plot(trialCount(csMinusTrials), smooth(TE.csLicks.rate(csMinusTrials), smoothFactor), 'r', 'LineWidth', 2);
% subplot(4,1,1); scatter(1:length(allTE.epoch), cuedReward.anticipatoryLicks2);
maxLR = max(TE.csLicks.rate(csPlusTrials));
% subplot(4,1,1); plot(smooth(cuedReward.anticipatoryLicks1 / mean(cuedReward.blLicks)));
% maxLR = smooth(max(cuedReward.anticipatoryLicks2 / mean(cuedReward.blLicks)));
stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxLR, 'g', 'Marker', 'none');
stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxLR, 'r', 'Marker', 'none');
ylabel('antic. licks'); title('CS+ trials');
set(gca, 'XLim', [1 length(trialCount)]);

subplot(2,1,2); hold on
scatter(trialCount(csMinusTrials), TE.csLicks.rate(csMinusTrials), 10, 'o'); % both blLicks and anticipatoryLicks2 are from 2s periods
plot(trialCount(csMinusTrials), smooth(TE.csLicks.rate(csMinusTrials), smoothFactor), 'b', 'LineWidth', 2);
% subplot(4,1,1); scatter(1:length(allTE.epoch), cuedReward.anticipatoryLicks2);
maxLR = max(TE.csLicks.rate(csMinusTrials));
% subplot(4,1,1); plot(smooth(cuedReward.anticipatoryLicks1 / mean(cuedReward.blLicks)));
% maxLR = smooth(max(cuedReward.anticipatoryLicks2 / mean(cuedReward.blLicks)));
hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxLR, 'g', 'Marker', 'none');
hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxLR, 'r', 'Marker', 'none');
ylabel('antic. licks'); title('CS- trials');
set(gca, 'XLim', [1 length(trialCount)]);

if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end

%% auROC plots
saveName = [subjectName '_auROC'];
h=ensureFigure(saveName, 1); 
plot(trialCount, TE.auROC, 'b', 'LineWidth', 1);
hold on;
yLims = [min(TE.auROC); max(TE.auROC)];
line(repmat(find(TE.BlockChange)', 2, 1), repmat(yLims, 1, sum(TE.BlockChange)), 'Color', 'g');
line(repmat(find(TE.sessionChange)', 2, 1), repmat(yLims, 1, sum(TE.sessionChange)), 'Color', 'r');

%% test thing for sanity check
% sanity check
TE.randVar = randn(nTrials, 1);
TE.auROCnoise = calc_auROC_reversalsFromTE(TE, TE.randVar, csPlusTrials, csMinusTrials, 20, 1);
saveName = [subjectName '_auROCnoise'];
h=ensureFigure(saveName, 1); 
plot(trialCount, TE.auROCnoise, 'b', 'LineWidth', 1);
hold on;
yLims = [min(TE.auROCnoise); max(TE.auROCnoise)];
line(repmat(find(TE.BlockChange)', 2, 1), repmat(yLims, 1, sum(TE.BlockChange)), 'Color', 'g');
line(repmat(find(TE.sessionChange)', 2, 1), repmat(yLims, 1, sum(TE.sessionChange)), 'Color', 'r');

    %% post-reversal analysis, assuming you always have photometry in channel 1, but not necessarily in channel 2
    
    dataToPull = {...
        'csLicks', TE.csLicks.rate,...
        'usLicks', TE.usLicks.rate,...
        'auROC', TE.auROC,...
        'trialOutcome', TE.trialOutcome,...
        'trialType', TE.trialType,...
        'trialNumber', TE.trialNumber,...
        'filename', TE.filename,...
        'ReinforcementOutcome', TE.ReinforcementOutcome,...
        'OdorValveIndex', TE.OdorValveIndex,...
        };
    
    RE = struct();
    RE.csPlus = extractReversalsFromTE(TE, csPlusTrials, dataToPull);%, 'maxReversals', 1);
    RE.csMinus = extractReversalsFromTE(TE, csMinusTrials, dataToPull);%, 'maxReversals', 1);
    RE.csPlusReward = extractReversalsFromTE(TE, csPlusTrials & rewardTrials, dataToPull);%, 'maxReversals', 1);
    RE.allTrials = extractReversalsFromTE(TE, validTrials, dataToPull);%, 'maxReversals', 1);
    nReversals = size(RE.csPlus.csLicks.after, 1);
    
if saveOn
    save(fullfile(savepath, ['RE_' subjectName '.mat']), 'RE');
    disp(['*** saving: ' fullfile(savepath, ['RE_' subjectName '.mat']) ' ***']);
end

%%
ensureFigure('test', 1);

subplot(2,1,1); hold on
xData = [RE.csMinus.trialsBefore RE.csPlus.trialsAfter];
revCsLicks = [RE.csMinus.csLicks.before RE.csPlus.csLicks.after];
plot(xData, nanfastsmooth(nanmean(revCsLicks), 5,1));

xData_newCsMinus = [RE.csPlus.trialsBefore RE.csMinus.trialsAfter];
revCsLicks_newCsMinus = [RE.csPlus.csLicks.before RE.csMinus.csLicks.after];
plot(xData_newCsMinus, nanfastsmooth(nanmean(revCsLicks_newCsMinus), 5,1));



