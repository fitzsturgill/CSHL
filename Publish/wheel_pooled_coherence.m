%{ 
Cross spectrogram for each trial, align magnitude and phase relative to
time from last reward. Key here is to use a sufficiently long moving window
to selectively detect low frequency (~0.2-5Hz) coherence between ChAT GCAMP6s and 
DAT RCaMP1a.
 
Goals: 
1) Coherence w.r.t. reward for each animal, bar graphs
2) Coherence (not time-resolved) for each animal, plot
3) cross correlation
4) Radial plot with coherence phase and magnitude?
5) Summary bar graph across n = 6 mice, coherence per-reward vs post-reward

%}

DB = dbLoadExperiment('wheel');
saveOn = 1;
savepath = fullfile(DB.path, filesep, 'pooled', filesep);

na = length(DB.animals);

Photometry = 'PhotometryHF';
    
movingWin = [5 0.5]; % 10 0.5
tapers = [3 5]; % 5 9
fpass = [1/120 10];
fpass_avg = [0.2 3]; 
maxTime = 40;
whiten = true; % works either way, just removes 1/f^a trend, and coherence is at bottom of frequency band


%%
% to compile data:
cohData = struct(...
    'cxcg', [],...
    'cxcgS', [],...
    'coh', [],...
    'cohS', [],...
    'C', [],...
    'CS', [],...
    'tfr', []... % time from reward
    );
cohData = repmat(cohData, na, 1);

% first calculate coherence, and time from reward, store in structure
for counter = 1:na
    animal = DB.animals{counter};
    dbLoadAnimal(DB, animal);        
    Fs = TE.(Photometry).sampleRate;
    r1 = TE.(Photometry).data(1).ZS;
    r2 = TE.(Photometry).data(2).ZS;

    nTrials = size(r1, 1);
    cxcg = bpCalcCrossCoherence(r1', r2', Fs, 'trialave', false, 'tapers', tapers, 'fpass', fpass, 'movingwin', movingWin, 'whiten', whiten);    
    cxcg_shift = bpCalcCrossCoherence(circshift(r1, 1, 1)', r2', Fs, 'trialave', false, 'tapers', tapers, 'fpass', fpass, 'movingwin', movingWin, 'whiten', whiten);
    
    % also calculate coherence, not time-resolved    
    coh = bpCalcCoherence(r1', r2', Fs, 'trialave', true, 'tapers', tapers, 'fpass', fpass, 'whiten', whiten);
    coh.Cerr = max(coh.Cerr, 0); % you can't have negative coherence
    coh_shift = bpCalcCoherence(circshift(r1, 1, 1)', r2', Fs, 'trialave', true, 'tapers', tapers, 'fpass', fpass, 'whiten', whiten);
    coh_shift.Cerr = max(coh_shift.Cerr, 0);
    C = permute(cxcg.C, [2 1 3]);
    C_shift = permute(cxcg_shift.C, [2 1 3]);
    C = C(:,:);
    C_shift = C_shift(:,:);
    nPoints = length(cxcg.t);    
    
    if counter == 1
        [cxcg_all, cxcg_shift_all] = deal(NaN(size(cxcg.C, 2), size(cxcg.C, 1), length(DB.animals)));
    end
    timeMatrix = bpCalcTimeFromEvent(TE, 'Reward', 'trialStart', TE.TrialStartTimestamp, 'dataStart', repmat(cxcg.t(1), nTrials, 1), 'Fs', 1/movingWin(2), 'duration', range(cxcg.t) + movingWin(2));
    timeMatrix = min(timeMatrix, maxTime);
    timeMatrix = timeMatrix';
    timeMatrix = timeMatrix(:); % concatenate trials
    sz = size(C);
    
    cohData(counter).cxcg = cxcg;
    cohData(counter).cxcgS = cxcg_shift;
    cohData(counter).coh = coh;
    cohData(counter).cohS = coh_shift;    
    cohData(counter).C = C;
    cohData(counter).CS = C_shift;
    cohData(counter).tfr = timeMatrix;
end

save(fullfile(savepath, 'coherence_pooled.mat'), 'cohData');
disp(['*** saved: ' fullfile(savepath, 'coherence_pooled.mat') ' ***']);
    
%%
% 1) Coherence w.r.t. reward for each animal, bar graphs
% 2) Coherence (not time-resolved) for each animal, plot
% 3) cross correlation
% 4) Radial plot with coherence phase and magnitude?
% 5) Summary bar graph across n = 6 mice, coherence per-reward vs post-reward

load(fullfile(savepath, 'coherence_pooled.mat'), 'cohData');
disp(['*** loaded: ' fullfile(savepath, 'coherence_pooled.mat') ' ***']);

%%
fh = zeros(na,1);
for counter = 1:na    
    animal = DB.animals{counter};
    fh(counter) = ensureFigure(sprintf('wheel_coh_%s', animal), 1);
    mcPortraitFigSetup(gcf);
    % bin means
    if fpass_avg(1)
        ix1 = crossing(cohData(counter).cxcg.f - fpass_avg(1));
    else
        ix1 = 1; % in case you are starting from frequency 0
    end
    ix2 = crossing(cohData(counter).cxcg.f - fpass_avg(2));
    bins = [0:movingWin(2):maxTime];
    [xC_Means, xC_Errors, xC_timeFromReward] = binnedMeansXY(cohData(counter).tfr, mean(cohData(counter).C(ix1:ix2,:))', bins);
    [xC_Means_shift, xC_Errors_shift, ~] = binnedMeansXY(cohData(counter).tfr, mean(cohData(counter).CS(ix1:ix2,:))', bins);
    
    subplot(3,2,1); hold on; title('cross spectrogram aligned');
    textBox(animal);
    errorbar(xC_timeFromReward, xC_Means, xC_Errors, 'b');
    errorbar(xC_timeFromReward, xC_Means_shift, xC_Errors_shift, 'k');    
    set(gca, 'YLim', [0 1]);
    xlabel('time from reward(s)'); ylabel('Coherence');
    
    subplot(3,2,2); hold on; title('cross spectrogram collapsed');
    plot(cohData(counter).cxcg.f, mean(cohData(counter).C, 2), 'b');
    plot(cohData(counter).cxcg.f, mean(cohData(counter).CS, 2), 'k');
    set(gca, 'YLim', [0 1]);
    addStimulusPatch(gca, fpass_avg);
    xlabel('frequency'); ylabel('Coherence');
    
    subplot(3,2,3); hold on; title('Coherence');
    boundedline(cohData(counter).cohS.f, cohData(counter).cohS.C, abs(cohData(counter).cohS.Cerr' - cohData(counter).cohS.C), 'k');
    boundedline(cohData(counter).coh.f, cohData(counter).coh.C, abs(cohData(counter).coh.Cerr' - cohData(counter).coh.C), 'b');
    set(gca, 'YLim', [0 1]);
    xlabel('frequency'); ylabel('Coherence');
    
    subplot(3,2,4); hold on; title('Coherence');
    boundedline(cohData(counter).cohS.f, cohData(counter).cohS.C, abs(cohData(counter).cohS.Cerr' - cohData(counter).cohS.C), 'k');
    boundedline(cohData(counter).coh.f, cohData(counter).coh.C, abs(cohData(counter).coh.Cerr' - cohData(counter).coh.C), 'b');
    set(gca, 'XScale', 'log');
    set(gca, 'YLim', [0 1]);
    xlabel('frequency'); ylabel('Coherence');                   
end

%% write pdfs
h = waitbar(0, 'slowly writing pdfs');
pdfn = fullfile(savepath, 'wheel_coherence_summary.pdf');
for counter = 1:na
        if counter == 1
        export_fig(fh(counter),pdfn);  % write to pdf
    else
        export_fig(fh(counter),'-append',pdfn);  % write to pdf
    end
    waitbar(counter/na);
end
close(h);

%% bar graph for coherence, peri-reward, post-reward

figSize = [1 1];
bins = [0 5; 5 120]; % 0-5s = peri,   5-120s = post
[periReward, postReward, periShift, postShift] = deal(zeros(na,1));
for counter = 1:na    
    % bin means
    if fpass_avg(1)
        ix1 = crossing(cohData(counter).cxcg.f - fpass_avg(1));
    else
        ix1 = 1; % in case you are starting from frequency 0
    end
    ix2 = crossing(cohData(counter).cxcg.f - fpass_avg(2));

    % peri
    ixtime = (cohData(counter).tfr > bins(1,1)) & (cohData(counter).tfr <= bins(1,2));
    periReward(counter) = mean(mean(cohData(counter).C(ix1:ix2,ixtime)));
    periShift(counter) = mean(mean(cohData(counter).CS(ix1:ix2,ixtime)));
    % post
    ixtime = (cohData(counter).tfr > bins(2,1)) & (cohData(counter).tfr <= bins(2,2));
    postReward(counter) = mean(mean(cohData(counter).C(ix1:ix2,ixtime)));
    postShift(counter) = mean(mean(cohData(counter).CS(ix1:ix2,ixtime)));

%     [xC_Means, xC_Errors, xC_timeFromReward] = binnedMeansXY(cohData(counter).tfr, mean(cohData(counter).C(ix1:ix2,:))', bins);
%     [xC_Means_shift, xC_Errors_shift, ~] = binnedMeansXY(cohData(counter).tfr, mean(cohData(counter).CS(ix1:ix2,:))', bins);    
end

comp = {'peri_vs_post'; 'peri_vs_shift'; 'post_vs_shift'};
p = zeros(3,1);
n = zeros(3,1);
coh_stats = table(comp, p, n);
[~, coh_stats.p(1)] = ttest(periReward, postReward);
[~, coh_stats.p(2)] = ttest(periReward, periShift);
[~, coh_stats.p(3)] = ttest(postReward, postShift);
coh_stats.n(:) = length(periReward);



%%
saveName = 'wheel_coh_barGraph';
ensureFigure(saveName, 1);
axes; hold on;
data = [periShift postShift];
errorbar(1:2, mean(data), std(data) / size(data, 1), 'k');
data = [periReward postReward];
errorbar(1:2, mean(data), std(data) / size(data, 1), 'b');
set(gca, 'Xlim', [0.5 2.5], 'XTick', [1 2], 'XTickLabel', {'peri', 'post'});
ylabel('Coherence'); xlabel('Reward');

formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
end
       
    %% spectrograms and cross spectrogram, comparison with trial data
    movingwin = [5 0.5]; % 1 .1
    tapers = [3 5];     %3 5
    whiten = false;
    fpass = [0 10];
    Photometry = 'PhotometryHF';
    climfactor = 4;
    
    Fs = TE.(Photometry).sampleRate;        
    
    d1 = TE.(Photometry).data(1).ZS';
    d2 = TE.(Photometry).data(2).ZS';
    d1 = zscore(d1, 0, 1);
    d2 = zscore(d2, 0, 1);
    xdata = TE.(Photometry).xData;
        
    
    sg1 = bpCalcSpectrogram(d1, Fs, 'movingwin', movingwin, 'whiten', whiten, 'tapers', tapers, 'trialave', 0, 'fpass', fpass);
    sg2 = bpCalcSpectrogram(d2, Fs, 'movingwin', movingwin, 'whiten', whiten, 'tapers', tapers, 'trialave', 0, 'fpass', fpass);
    cxcg = bpCalcCrossCoherence(d1, d2, Fs, 'whiten', whiten, 'trialave', false, 'tapers', tapers, 'fpass', fpass, 'movingwin', movingwin);
    
    %%
    trial = 1;
    if whiten
        saveName = 'coherence_dev_whitened';
    else
        saveName = 'coherence_dev';
    end
    fig = ensureFigure(saveName, 1);
    mcLandscapeFigSetup(fig);
    fig.Color = [0.7 0.7 0.7];
    
    rewards = TE.Reward{trial} - TE.(Photometry).startTime(trial);
    rewards = rewards(:,1);    
    
    subplot(5,1,1); hold on;
    line(repmat(rewards', 2, 1), repmat([-3; 3], 1, length(rewards)), 'Color', 'b', 'LineWidth', 2);
    plot(xdata, d1(:,trial), 'g');
    plot(xdata, d2(:,trial), 'r');

    addOrginLines;
    ylabel('Fluor (Zscore)');
    
    subplot(5,1,2); hold on;
    plot(xdata(1:end-1), diff(d1(:,trial)), 'g');
    plot(xdata(1:end-1), diff(d2(:,trial)), 'r');
    set(gca, 'YLim', [-0.1 0.1]);
    ylabel('Derivative');
    
    
    ax2 = subplot(5,1,3);
    cdata = squeeze(sg1.S(:,:,trial))';
    clim = [mean(cdata(:)) - std(cdata(:)*climfactor) mean(cdata(:)) + std(cdata(:))*climfactor];
    imagesc([sg1.t(1) sg1.t(end)], [sg1.f(1) sg1.f(end)], cdata, clim); set(gca, 'YDir', 'normal');
    ax2.YColor = [0 1 0];
    c = colorbar; c.Label.String = 'Ch1. Power'; c.Location = 'north'; 
    c.Color = [0.9 0.9 0.9];    
    ylabel('Frequency');
    
    ax3 = subplot(5,1,4); 
    cdata = squeeze(sg2.S(:,:,trial))';
    clim = [mean(cdata(:)) - std(cdata(:))*climfactor mean(cdata(:)) + std(cdata(:))*climfactor];
    imagesc([sg1.t(1) sg1.t(end)], [sg1.f(1) sg1.f(end)], cdata, clim); set(gca, 'YDir', 'normal'); 
    ax3.YColor = [1 0 0];
    c = colorbar; c.Label.String = 'Ch2. Power'; c.Location = 'north'; 
    c.Color = [0.9 0.9 0.9];
    ylabel('Frequency');
        
    subplot(5,1,5);
    cdata = squeeze(cxcg.C(:,:,trial))';
    clim = [mean(cdata(:)) - std(cdata(:))*climfactor mean(cdata(:)) + std(cdata(:))*climfactor];
    imagesc([cxcg.t(1) cxcg.t(end)],[cxcg.f(1) cxcg.f(end)], cdata, clim);
    set(gca, 'YDir', 'normal');
    c = colorbar; c.Label.String = 'Magnitude Squared Coherence'; c.Location = 'north';
    c.Color = [0.9 0.9 0.9];
    xlabel('Time (s)'); ylabel('Frequency');
    
    if saveOn    
        print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));
    end
    
    

