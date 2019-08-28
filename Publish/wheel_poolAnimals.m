% wheel_poolAnimals

DB = dbLoadExperiment('wheel');
saveOn = 1;
savepath = fullfile(DB.path, filesep, 'pooled', filesep);

na = length(DB.animals);
%%
%{ 
goals:
1) noise correlations for reward responses
2) coherence
%}


window = [0.1 0.6];
bl_window = [-0.5 0];
maxLag = 5;
Rnoise = struct(...
    'rew', NaN(na, 1),...
    'rew_shift', NaN(na, 1),...
    'rew_delta', NaN(na, 1),...
    'rew_delta_shift', NaN(na, 1)...
    );

for counter = 1:na
    animal = DB.animals{counter};
    dbLoadAnimal(DB, animal);
        
    r1 = extractDataByTimeStamps(TE.Photometry.data(1).ZS, TE.Photometry.startTime, 20, TE.Reward, window);
    if strcmp(animal, 'DC_36') % kludge to deal with weird delayed reward response for DC_36
        r2 = extractDataByTimeStamps(TE.Photometry.data(2).ZS, TE.Photometry.startTime, 20, TE.Reward, window + 0.5); 
    else
        r2 = extractDataByTimeStamps(TE.Photometry.data(2).ZS, TE.Photometry.startTime, 20, TE.Reward, window);
    end
    b1 = extractDataByTimeStamps(TE.Photometry.data(1).ZS, TE.Photometry.startTime, 20, TE.Reward, bl_window);    
    b2 = extractDataByTimeStamps(TE.Photometry.data(2).ZS, TE.Photometry.startTime, 20, TE.Reward, bl_window);
    
    goodones = all(isfinite(r1), 2) & all(isfinite(r2), 2) & all(isfinite(b1), 2) & all(isfinite(b2), 2);
    r1 = r1(goodones,:);
    r2 = r2(goodones,:);
    b1 = b1(goodones,:);
    b2 = b2(goodones,:);
    
    r2_noise = nanmean(r2, 2);
    r2_noise = r2_noise - nanmean(r2_noise);
    r1_noise = nanmean(r1, 2);
    r1_noise = r1_noise - nanmean(r1_noise);
    
    d1 = nanmean(r1,2) - nanmean(b1,2);
    d2 = nanmean(r2,2) - nanmean(b2,2);
    
    Rnoise.rew(counter) = corr(r1_noise, r2_noise);
    Rnoise.rew_shift(counter) = corr(r1_noise, circshift(r2_noise, 1));
    Rnoise.rew_delta(counter) = corr(d1, d2);
    Rnoise.rew_delta_shift(counter) = corr(d1, circshift(d2, 1));   
end

%%
Photometry = 'PhotometryHF';
xcWindow = [-2 10];
movingWin = [5 0.5];
tapers = [3 5];
fpass = [0 10];
saveName = 'coherence';
ensureFigure(saveName, 1); 
mcLandscapeFigSetup(gcf);
for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    dbLoadAnimal(DB, animal);
        
    r1 = extractDataByTimeStamps(TE.(Photometry).data(1).ZS, TE.(Photometry).startTime, TE.(Photometry).sampleRate, TE.Reward, xcWindow);
    r2 = extractDataByTimeStamps(TE.(Photometry).data(2).ZS, TE.(Photometry).startTime, TE.(Photometry).sampleRate, TE.Reward, xcWindow);
%     r1 = diff(r1, 1, 2);
%     r2 = diff(r2, 1, 2);
    cxcg = bpCalcCrossCoherence(r1', r2', TE.(Photometry).sampleRate, 'trialave', false, 'tapers', tapers, 'fpass', fpass, 'movingwin', movingWin, 'whiten', false);
    cxcg_shift = bpCalcCrossCoherence(circshift(r1, 1, 1)', r2', TE.(Photometry).sampleRate, 'trialave', false, 'tapers', tapers, 'fpass', fpass, 'movingwin', movingWin, 'whiten', false);
    if counter == 1
        [cxcg_all, cxcg_shift_all, phase_all, phase_shift_all] = deal(NaN(size(cxcg.C, 2), size(cxcg.C, 1), length(DB.animals)));
    end
    cxcg_all(:,:,counter) = squeeze(nanmean(cxcg.C, 3))';
    cxcg_shift_all(:,:,counter) = squeeze(nanmean(cxcg_shift.C, 3))';
    phase_all(:,:,counter) = squeeze(nanmean(cxcg.phi, 3))';
    phase_shift_all(:,:,counter) = squeeze(nanmean(cxcg_shift.phi, 3))';
    subplot(2,3,counter); title(animal);
    imagesc([cxcg.t(1) cxcg.t(end)],[cxcg.f(1) cxcg.f(end)], squeeze(nanmean(cxcg.C, 3))', [0 1]);
    textBox(animal);
    set(gca, 'YDir', 'normal');
    ylabel('Frequency'); xlabel('Time from rew. (s)');   
end

if saveOn    
%     print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));
end

%% 
ensureFigure('test4', 1);
mn = 2;
subplot(2,2,1); imagesc(squeeze(cxcg_all(:,:,mn)));
subplot(2,2,2); imagesc(squeeze(phase_all(:,:,mn)));
subplot(2,2,3); imagesc(squeeze(cxcg_shift_all(:,:,mn)));
subplot(2,2,4); imagesc(squeeze(phase_shift_all(:,:,mn)));
%%
figSize = [2 1];
clim = [0.3 1];
saveName = 'coherencegram_pooled';
ensureFigure(saveName, 1);
imagesc([cxcg.t(1) cxcg.t(end)],[cxcg.f(1) cxcg.f(end)], squeeze(nanmean(cxcg_all, 3)), clim);
xlabel('time from rew. (s)');
ylabel('F (hz)');
set(gca, 'YDir', 'normal'); 
set(gca, 'YLim', [0 6]);
colorbar;
formatFigurePublish('size', figSize);
if saveOn    
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));    
end

saveName = 'coherencegram_pooled_shifted';
ensureFigure(saveName, 1);
imagesc([cxcg_shift.t(1) cxcg_shift.t(end)],[cxcg_shift.f(1) cxcg_shift.f(end)], squeeze(nanmean(cxcg_shift_all, 3)), clim);
xlabel('time from rew. (s)');
ylabel('F (hz)');
set(gca, 'YDir', 'normal'); 
set(gca, 'YLim', [0 6]);
colorbar;
formatFigurePublish('size', figSize);
if saveOn    
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));    
end


%% plot rasters for each mouse and save them together in a pdf
climfactor = 2;
fh = zeros(na,1);
h = waitbar(0, 'slowly writing pdfs');
pdfn = fullfile(DB.path, 'pooled', 'wheel_rasters.pdf');
for counter = 1:na
    animal = DB.animals{counter};
    dbLoadAnimal(DB, animal);
    
    window = [-2 2];
    [rewards_dat, ts, tn] = extractDataByTimeStamps(TE.Photometry.data(2).ZS, TE.Photometry.startTime, 20, TE.Reward, [-2 2]);
    rewards_chat = extractDataByTimeStamps(TE.Photometry.data(1).ZS, TE.Photometry.startTime, 20, TE.Reward, [-2 2]);
    rewards_pupil = extractDataByTimeStamps(TE.pupil.pupDiameterNorm, TE.Photometry.startTime, 20, TE.Reward, [-2 2]);
    fh(counter) = ensureFigure(sprintf('rasters_%s', animal), 1);     
    mcPortraitFigSetup(fh(counter));
    subplot(3,2,1); triggeredEventRasterFromTE(TE, 'Port1In', TE.Reward); textBox(animal, [], [], 16);
    subplot(3,2,2); image(rewards_pupil, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', [nanmean(rewards_pupil(:)) - nanstd(rewards_pupil(:)) * climfactor nanmean(rewards_pupil(:)) + nanstd(rewards_pupil(:)) * climfactor]); colormap('parula');  title('Pupil');    
    subplot(3,2,3); image(rewards_chat, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', [nanmean(rewards_chat(:)) - nanstd(rewards_chat(:)) * climfactor nanmean(rewards_chat(:)) + nanstd(rewards_chat(:)) * climfactor]); colormap('parula');  title('ChAT');
    subplot(3,2,4); image(rewards_dat, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', [nanmean(rewards_dat(:)) - nanstd(rewards_dat(:)) * climfactor nanmean(rewards_dat(:)) + nanstd(rewards_dat(:)) * climfactor]); colormap('parula');    title('DAT');
    xdata = linspace(window(1), window(2), size(rewards_chat, 2));
    subplot(3,2,5); plot(xdata, nanmean(rewards_chat));
    xlabel('time from reward (s)');
    subplot(3,2,6); plot(xdata, nanmean(rewards_dat));
    xlabel('time from reward (s)');
    
    if counter == 1
        export_fig(gcf,pdfn);  % write to pdf
    else
        export_fig(gcf,'-append',pdfn);  % write to pdf
    end
    waitbar(counter/na);
end
close(h);





