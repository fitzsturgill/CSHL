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

%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end
%%  generate trial lookups for different combinations of conditions
LNL_conditions;

%% calculate discriminability of csPlus from csMinus trials according to anticipatory licking
[D, P, CI] = calc_auROC_reversalsFromTE(TE, TE.csLicks.rate, csPlusTrials, csMinusTrials,...
    'window', 30, 'windowMode', 'global', 'reset', 1, 'nBoot', 1000);
%
% [data1, data2, ~] = calc_auROC_reversalsFromTE(TE, TE.csLicks.rate, csPlusTrials, csMinusTrials,...
%     'window', 10, 'reset', 1, 'nBoot', 100);

TE.auROC = struct('D', [], 'P', [], 'CI', []);
TE.auROC.D = D;
TE.auROC.P = P;
TE.auROC.CI = CI;



%% auROC plots
saveName = [subjectName '_auROCboot'];
h=ensureFigure(saveName, 1); 
% plot(trialCount, TE.auROC.D, '-k', 'LineWidth', 1);
hold on;
colormap jet;
scatter(trialCount, TE.auROC.D, 20, TE.auROC.P);
plot(trialCount, TE.auROC.P);
yLims1 = [min(TE.auROC.D); max(TE.auROC.D)];
yLims2 = [min(TE.auROC.D) * 0.9; max(TE.auROC.D) * 0.9];
line(repmat(find(TE.sessionChange)', 2, 1), repmat(yLims1, 1, sum(TE.sessionChange)), 'Color', 'r');
line(repmat(find(TE.BlockChange)', 2, 1), repmat(yLims2, 1, sum(TE.BlockChange)), 'Color', 'g');



%% auROC plots confidence interval
saveName = [subjectName '_auROCboot_CI2'];
h=ensureFigure(saveName, 1); 

boundedline(trialCount(:), nan2zero(TE.auROC.D(:)), abs([nan2val(TE.auROC.CI(:,1), 0) nan2val(TE.auROC.CI(:,2), 0)] - nan2zero(TE.auROC.D(:))), 'k', 'alpha', 'transparency', 0.1); hold on;
hold on; colormap jet
scatter(trialCount, TE.auROC.D, 20, TE.auROC.P);


yLims1 = [min(TE.auROC.D); max(TE.auROC.D)];
yLims2 = [min(TE.auROC.D) * 0.9; max(TE.auROC.D) * 0.9];
line(repmat(find(TE.sessionChange)', 2, 1), repmat(yLims1, 1, sum(TE.sessionChange)), 'Color', 'r');
line(repmat(find(TE.BlockChange)', 2, 1), repmat(yLims2, 1, sum(TE.BlockChange)), 'Color', 'g');
addOrginLines;
set(gca, 'YLim', [-1.1 1.1]);

%% plot p value
saveName = [subjectName '_auROCboot_pValue_global'];
h=ensureFigure(saveName, 1);
scatter(trialCount(:), TE.auROC.P, 20, TE.auROC.D);
yLims1 = [0 1];
yLims2 = [0.1; 0.9];
% line(repmat(find(TE.sessionChange)', 2, 1), repmat(yLims1, 1, sum(TE.sessionChange)), 'Color', 'r');
line(repmat(find(TE.BlockChange)', 2, 1), repmat(yLims2, 1, sum(TE.BlockChange)), 'Color', 'g');
addOrginLines;


%% test thing for sanity check
% sanity check
TE.randVar = randn(nTrials, 1);
[D, P, CI] = calc_auROC_reversalsFromTE(TE, TE.randVar, csPlusTrials, csMinusTrials,...
    'window', 20, 'reset', 1, 'nBoot', 100);
TE.auROCnoise = [];
TE.auROCnoise.D = D;
TE.auROCnoise.P = P;
TE.auROCnoise.CI = CI;
saveName = [subjectName '_auROCnoise'];
h=ensureFigure(saveName, 1); 

% boundedline(trialCount(:), nan2zero(TE.auROCnoise.D(:)), nan2zero(TE.auROCnoise.P(:)), 'b'); hold on;
% hold on;
% yLims1 = [min(TE.auROCnoise.D); max(TE.auROCnoise.D)];
% yLims2 = [min(TE.auROCnoise.D) * 0.9; max(TE.auROCnoise.D) * 0.9];
% line(repmat(find(TE.sessionChange)', 2, 1), repmat(yLims, 1, sum(TE.sessionChange)), 'Color', 'r');
% line(repmat(find(TE.BlockChange)', 2, 1), repmat(yLims, 1, sum(TE.BlockChange)), 'Color', 'g');
% addOrginLines;

boundedline(trialCount(:), nan2zero(TE.auROCnoise.D(:)), abs([nan2val(TE.auROCnoise.CI(:,1), 0) nan2val(TE.auROCnoise.CI(:,2), 0)] - nan2zero(TE.auROCnoise.D(:))), 'k', 'alpha', 'transparency', 0.1); hold on;
hold on; colormap jet
scatter(trialCount, TE.auROCnoise.D, 20, TE.auROCnoise.P);


yLims1 = [min(TE.auROCnoise.D); max(TE.auROCnoise.D)];
yLims2 = [min(TE.auROCnoise.D) * 0.9; max(TE.auROCnoise.D) * 0.9];
line(repmat(find(TE.sessionChange)', 2, 1), repmat(yLims1, 1, sum(TE.sessionChange)), 'Color', 'r');
line(repmat(find(TE.BlockChange)', 2, 1), repmat(yLims2, 1, sum(TE.BlockChange)), 'Color', 'g');
addOrginLines;
set(gca, 'YLim', [-1.1 1.1]);



    %% post-reversal analysis, assuming you always have photometry in channel 1, but not necessarily in channel 2
    
    dataToPull = {...
        'csLicks', TE.csLicks.rate,...
        'usLicks', TE.usLicks.rate,...
        'auROC_D', TE.auROC.D,...
        'auROC_P', TE.auROC.P,...        
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

%%
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

