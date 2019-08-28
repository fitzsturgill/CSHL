%% different strategy, cross spectrogram for all the data, don't average then unwrap relative to time from most recent reward

DB = dbLoadExperiment('wheel');
saveOn = 1;
savepath = fullfile(DB.path, filesep, 'pooled', filesep);

na = length(DB.animals);

Photometry = 'PhotometryHF';
    animalnumber = 1;
    animal = DB.animals{animalnumber};

    
    
movingWin = [5 0.5]; % 10 0.5
tapers = [3 5]; % 5 9
fpass = [0 10];
fpass_avg = [0.2 3]; 
maxTime = 40;

whiten = true;


ensureFigure('binned_coherence', 1);
for counter = 1:na
    animal = DB.animals{counter};
    dbLoadAnimal(DB, animal);        
    Fs = TE.(Photometry).sampleRate;
    r1 = TE.(Photometry).data(1).ZS;    
    r2 = TE.(Photometry).data(2).ZS;

    nTrials = size(r1, 1);
    cxcg = bpCalcCrossCoherence(r1', r2', Fs, 'trialave', false, 'tapers', tapers, 'fpass', fpass, 'movingwin', movingWin, 'whiten', whiten);
    cxcg_shift = bpCalcCrossCoherence(circshift(r1, 1, 1)', r2', Fs, 'trialave', false, 'tapers', tapers, 'fpass', fpass, 'movingwin', movingWin, 'whiten', whiten);
    %
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

    
    
    
    % bin means
    if fpass_avg(1)
        ix1 = crossing(cxcg.f - fpass_avg(1));
    else
        ix1 = 1; % in case you are starting from frequency 0
    end
    ix2 = crossing(cxcg.f - fpass_avg(2));
    bins = [0:movingWin(2):maxTime];
    [xC_Means, xC_Errors, xC_timeFromReward] = binnedMeansXY(timeMatrix, mean(C(ix1:ix2,:))', bins);
    [xC_Means_shift, xC_Errors_shift, ~] = binnedMeansXY(timeMatrix, mean(C_shift(ix1:ix2,:))', bins);
    
    subplot(3,2,counter); hold on;
    textBox(animal);
    errorbar(xC_timeFromReward, xC_Means, xC_Errors, 'k');
    errorbar(xC_timeFromReward, xC_Means_shift, xC_Errors_shift, 'Color', [0.7 0.7 0.7]);
    xlabel('time from reward(s)'); ylabel('Coherence');
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
    
    

