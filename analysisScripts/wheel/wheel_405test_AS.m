saveOn = 0;
%%
saveOn = 1;
%%
sessions = bpLoadSessions; % load sessions
%% 
TE = makeTE_wheel_alternateLEDs(sessions); % make TE structure

%% Now saved in directory according to first session filename
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
basepath = 'Z:\SummaryAnalyses\wheel_405test\';
if length(unique(TE.filename)) > 1
    sep = strfind(TE.filename{1}, '_');
    subjectName = TE.filename{1}(1:sep(2)-1);
else
    sep = strfind(TE.filename{1}, '.');
    subjectName = TE.filename{1}(1:sep(1)-1);
end
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);


%%

dFFMode = {'simple', 'simple'};
try 
    bl = [0 sessions(1).SessionData.NidaqData{1,2}.duration(1) - 0.2];
catch
    bl = [0 9.9];
end

% channel 1 demodulated via channel 1 reference
TE.Photometry_ch1 = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', 1, 'baseline', bl, 'startField', 'Baseline', 'downsample', 305, 'forceAmp', 1);

% channel 1 demodulated via channel 2 reference
TE.Photometry_ch2 = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', 1, 'baseline', bl, 'startField', 'Baseline', 'downsample', 305, 'refChannels', 2, 'forceAmp', 1);


if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end
%%

LED1trials = TE.LED1_amp > 0 & TE.LED2_amp == 0;
LED2trials = TE.LED2_amp > 0 & TE.LED1_amp == 0;
LED12trials = TE.LED1_amp > 0 & TE.LED2_amp > 0;

%% do bleach fit manually on a trial by trial basis
nTrials = length(TE.filename);
nPoints = size(TE.Photometry_ch1.data(1).raw, 2);
TE.Photometry_ch1.data(1).bleach_fit = zeros(nTrials, nPoints);
TE.Photometry_ch2.data(1).bleach_fit = zeros(nTrials, nPoints);
TE.Photometry_ch1.data(1).bleach_dF = zeros(nTrials, nPoints);
TE.Photometry_ch2.data(1).bleach_dF = zeros(nTrials, nPoints);

x = (0:nPoints - 1);
startFitTime = 1; % start bleach fit n seconds into each trial
Fs = TE.Photometry_ch1.sampleRate;
for counter = 1:nTrials
    if TE.LED1_amp(counter)
        trialData = TE.Photometry_ch1.data(1).raw(counter, :);
        if any(isnan(trialData))
            trialData = inpaint_nans(trialData);
        end
        fo = fitoptions('Method', 'NonlinearLeastSquares',...
            'Upper', [Inf Inf 0],...
            'Lower', [0 -Inf -1],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
            'StartPoint', [min(trialData) range(trialData)/2 -5]...
            );
        model = 'a + b*exp(c*x)';
        ft = fittype(model, 'options', fo);        
        [fitobject, gof, output] = ... % new
            fit(x(bpX2pnt(startFitTime, Fs):end)', trialData(bpX2pnt(startFitTime, Fs):end)', ft, fo);        
        trialFit = fitobject.a + fitobject.b * exp(fitobject.c * x);
        TE.Photometry_ch1.data(1).bleach_fit(counter,:) = trialFit;
        TE.Photometry_ch1.data(1).bleach_dF(counter,:) = trialData - trialFit;
    end
    
    if TE.LED2_amp(counter)
        trialData = TE.Photometry_ch2.data(1).raw(counter, :);
        if any(isnan(trialData))
            trialData = inpaint_nans(trialData);
        end        
        fo = fitoptions('Method', 'NonlinearLeastSquares',...
            'Upper', [Inf Inf 0],...
            'Lower', [0 -Inf -1],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
            'StartPoint', [min(trialData) range(trialData)/2 -5]...
            );
        model = 'a + b*exp(c*x)';
        ft = fittype(model, 'options', fo);                
        [fitobject, gof, output] = ... % new
            fit(x(bpX2pnt(startFitTime, Fs):end)', trialData(bpX2pnt(startFitTime, Fs):end)', ft, fo);        
        trialFit = fitobject.a + fitobject.b * exp(fitobject.c * x);
        TE.Photometry_ch2.data(1).bleach_fit(counter,:) = trialFit;
        TE.Photometry_ch2.data(1).bleach_dF(counter,:) = trialData - trialFit;
    end    
end

%%
% nPoints = numel(TE.Photometry_ch1.data(1).raw);
% x = 0:nPoints-1;
% fo = TE.Photometry_ch1.bleachFit.fitobject_session;
% ch1_blFit = fo.a + fo.b * exp(fo.c * x);% + fo.d * exp(fo.e * x);
% fo = TE.Photometry_ch2.bleachFit.fitobject_session;
% ch2_blFit = fo.a + fo.b * exp(fo.c * x);% + fo.d * exp(fo.e * x);
% 
% saveName = 'expFitConcat_bleachFit';
% ensureFigure(saveName, 1);
% subplot(1,1,1); 
% ch1data = TE.Photometry_ch1.data(1).raw';
% ch1data = ch1data(:);
% plot(ch1data, 'b'); hold on;
% plot(ch1_blFit, 'c', 'LineWidth', 2);
% title('470nm');
% % subplot(1,2,2); 
% yyaxis right;
% ch2data = TE.Photometry_ch2.data(1).raw';
% ch2data = ch2data(:);
% plot(ch2data, 'm'); hold on;
% plot(ch2_blFit, 'r', 'LineWidth', 2);
% % title('405nm');
% legend('470', '470fit', '405', '405fit'); legend('boxoff');
% 
%     if saveOn
%         saveas(gcf, fullfile(savepath, saveName), 'fig');
%         saveas(gcf, fullfile(savepath, saveName), 'jpeg');
%     end
    
%%    

window = [-2 2];
fs = 20; % sample rate
blSamples = (0 - window(1)) * fs;

% LED 470 and 405
[rewards_ch1_both, ts, tn] = extractDataByTimeStamps(TE.Photometry_ch1.data(1).raw, TE.Photometry_ch1.startTime, 20, TE.Reward, [-2 2], LED12trials);
rewards_ch2_both = extractDataByTimeStamps(TE.Photometry_ch2.data(1).raw, TE.Photometry_ch1.startTime, 20, TE.Reward, [-2 2], LED12trials);
bl_ch1_both = nanmean(rewards_ch1_both(:,1:blSamples), 2);
rewards_ch1_both = bsxfun(@minus, rewards_ch1_both, bl_ch1_both);
bl_ch2_both = nanmean(rewards_ch2_both(:,1:blSamples), 2);
rewards_ch2_both = bsxfun(@minus, rewards_ch2_both, bl_ch2_both);

% LED 470 only
[rewards_ch1_470, ts, tn] = extractDataByTimeStamps(TE.Photometry_ch1.data(1).raw, TE.Photometry_ch2.startTime, 20, TE.Reward, [-2 2], LED1trials);
rewards_ch2_470 = extractDataByTimeStamps(TE.Photometry_ch2.data(1).raw, TE.Photometry_ch2.startTime, 20, TE.Reward, [-2 2], LED1trials);
bl_ch1_470 = nanmean(rewards_ch1_470(:,1:blSamples), 2);
rewards_ch1_470 = bsxfun(@minus, rewards_ch1_470, bl_ch1_470);
bl_ch2_470 = nanmean(rewards_ch2_470(:,1:blSamples), 2);
rewards_ch2_470 = bsxfun(@minus, rewards_ch2_470, bl_ch2_470);

% LED 405 only
[rewards_ch1_405, ts, tn] = extractDataByTimeStamps(TE.Photometry_ch1.data(1).raw, TE.Photometry_ch1.startTime, 20, TE.Reward, [-2 2], LED2trials);
rewards_ch2_405 = extractDataByTimeStamps(TE.Photometry_ch2.data(1).raw, TE.Photometry_ch2.startTime, 20, TE.Reward, [-2 2], LED2trials);
bl_ch1_405 = nanmean(rewards_ch1_405(:,1:blSamples), 2);
rewards_ch1_405 = bsxfun(@minus, rewards_ch1_405, bl_ch1_405);
bl_ch2_405 = nanmean(rewards_ch2_405(:,1:blSamples), 2);
rewards_ch2_405 = bsxfun(@minus, rewards_ch2_405, bl_ch2_405);


% calc mean and SD for ch1 and ch2 (both LEDs on), use this across all LED
% conditions
sd_ch1 = nanmean(nanstd(rewards_ch1_both(:,1:blSamples)));
sd_ch2 = nanmean(nanstd(rewards_ch2_both(:,1:blSamples)));

%% stitch together 
% nTrials = length(TE.filename);
% nSamples = size(rewards_dat_12, 2);
% rewards_dat_all = NaN(nTrials, nSamples);
% rewards_chat_all = NaN(nTrials, nSamples);
% 
% rewards_dat_all(LED12trials, :) = rewards_dat_12(:,:);
% rewards_dat_all(LED1trials, :) = rewards_dat_1(:,:);
% rewards_dat_all(LED2trials, :) = rewards_dat_2(:,:);
% 
% rewards_chat_all(LED12trials, :) = rewards_chat_12(:,:);
% rewards_chat_all(LED1trials, :) = rewards_chat_1(:,:);
% rewards_chat_all(LED2trials, :) = rewards_chat_2(:,:);



climFactor = 4;
clim_ch1 = [-climFactor * sd_ch1, climFactor * sd_ch1];
clim_ch2 = [-climFactor * sd_ch2, climFactor * sd_ch2];

ensureFigure('random_rewards_alternateLEDs', 1); 

% first column, ch1
subplot(3,2,1); image(rewards_ch1_both, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_ch1); colormap('jet');  ylabel('deltaF, both LEDs'); title('Ch 1');


subplot(3,2,3); image(rewards_ch1_470, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_ch1); colormap('jet');  ylabel('LED1 only');


subplot(3,2,5); image(rewards_ch1_405, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_ch1); colormap('jet');  ylabel('LED2 only');

% second column, ch2
subplot(3,2,2); image(rewards_ch2_both, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_ch2); colormap('jet');    title('Ch 2');


subplot(3,2,4); image(rewards_ch2_470, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_ch2); colormap('jet');    


subplot(3,2,6); image(rewards_ch2_405, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_ch2); colormap('jet');    



if saveOn
    saveas(gcf, fullfile(savepath, 'random_rewards_alternateLEDs.fig'));
    saveas(gcf, fullfile(savepath, 'random_rewards_alternateLEDs.jpg'));
end

%% save plot of baseline fluorecense values across conditions
fbl = zeros(1,6);
fbl_sem = zeros(1,6);
fbl(1) = nanmean(TE.Photometry_ch1.data(1).blF_raw(LED12trials));
fbl(2) = nanmean(TE.Photometry_ch1.data(1).blF_raw(LED1trials));
fbl(3) = nanmean(TE.Photometry_ch1.data(1).blF_raw(LED2trials));
fbl(4) = nanmean(TE.Photometry_ch2.data(1).blF_raw(LED12trials));
fbl(5) = nanmean(TE.Photometry_ch2.data(1).blF_raw(LED1trials));
fbl(6) = nanmean(TE.Photometry_ch2.data(1).blF_raw(LED2trials));

fbl_sem(1) = nanSEM(TE.Photometry_ch1.data(1).blF_raw(LED12trials));
fbl_sem(2) = nanSEM(TE.Photometry_ch1.data(1).blF_raw(LED1trials));
fbl_sem(3) = nanSEM(TE.Photometry_ch1.data(1).blF_raw(LED2trials));
fbl_sem(4) = nanSEM(TE.Photometry_ch2.data(1).blF_raw(LED12trials));
fbl_sem(5) = nanSEM(TE.Photometry_ch2.data(1).blF_raw(LED1trials));
fbl_sem(6) = nanSEM(TE.Photometry_ch2.data(1).blF_raw(LED2trials));

% fbl(1) = nanmean(TE.Photometry.data(1).blF_raw(LED12trials));
% fbl(2) = nanmean(TE.Photometry.data(2).blF_raw(LED12trials));
% fbl(3) = nanmean(TE.Photometry.data(1).blF_raw(LED1trials));
% fbl(4) = nanmean(TE.Photometry_1alt.data(2).blF_raw(LED1trials));
% fbl(5) = nanmean(TE.Photometry_2alt.data(1).blF_raw(LED2trials));
% fbl(6) = nanmean(TE.Photometry.data(2).blF_raw(LED2trials));
% 
% fbl_sem(1) = nanSEM(TE.Photometry.data(1).blF_raw(LED12trials));
% fbl_sem(2) = nanSEM(TE.Photometry.data(2).blF_raw(LED12trials));
% fbl_sem(3) = nanSEM(TE.Photometry.data(1).blF_raw(LED1trials));
% fbl_sem(4) = nanSEM(TE.Photometry_1alt.data(2).blF_raw(LED1trials));
% fbl_sem(5) = nanSEM(TE.Photometry_2alt.data(1).blF_raw(LED2trials));
% fbl_sem(6) = nanSEM(TE.Photometry.data(2).blF_raw(LED2trials));
ensureFigure('baseline_fluor', 1);

bar(1:6, fbl, 'b'); hold on;
errorbar(1:6,fbl, fbl_sem, '.', 'LineWidth', 2 );
set(gca, 'XTick', 1:6, 'XTickLabel',...
    {'ch1 both', 'ch1 470', 'ch1 405', 'ch2 both', 'ch2 470', 'ch2 405'});
ylabel('Baseline Fluor (V)');
if saveOn
    saveas(gcf, fullfile(savepath, 'baseline_fluor.fig'));
    saveas(gcf, fullfile(savepath, 'baseline_fluor.jpg'));
end

%% given that you've used 