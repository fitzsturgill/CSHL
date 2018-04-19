 %% script to run after LNL_odor_V2_pav_rev_AS
 
 
 %%
    ensureFigure('all_behavior_CsPlus', 1);
    reversals = find(diff(TE.BlockNumber(csPlusTrials, :))) + 1;
    subplot(1,4,1);
    
    image(TE.Wheel.data.V(csPlusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
    set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 2, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 2]); 
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('Velocity');
    subplot(1,4,2);
    image(TE.pupil.pupDiameterNorm(csPlusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
    set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    colormap('parula');  
    title('Pupil Diameter');    
    subplot(1,4,3); phRasterFromTE(TE, csPlusTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('ChAT'); xlabel('Time frome odor (s)');
    subplot(1,4,4); phRasterFromTE(TE, csPlusTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('DAT');    
    ax = findobj(gcf, 'Type', 'Axes');
%     set(ax, 'YLim', [40 80]);
    

%% allBehavior_rasters_SFN
    firstTrial = 100;
    all_behavior_trials = csPlusTrials & rewardTrials & hitTrials & ~isnan(mean(TE.pupil.pupDiameterNorm, 2)) & TE.trialNumber >= firstTrial;
    
    savename = 'allBehavior_rasters';
    ensureFigure(savename, 1);
%     reversals = find(diff(TE.BlockNumber(csPlusTrials, :))) + 1;
    
    subplot(1,5,1); 
    eventRasterFromTE(TE, all_behavior_trials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('Licks'); ylabel('trial number');
    set(gca, 'XLim', [-4 7]);
    
%     set(gca, 'YLim', [0 max(trialCount)]);
%     set(gca, 'FontSize', 14)
    wp = [bpX2pnt(-4, 20, -4) bpX2pnt(3, 20, -4)];
    subplot(1,5,2);    
    image(TE.Wheel.data.V(all_behavior_trials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
    set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 2, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 2]); 
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('Velocity');
    set(gca, 'YTick', []); 

    subplot(1,5,3);
    image(TE.pupil.pupDiameterNorm(all_behavior_trials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
    set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
    set(gca, 'XLim', [-4 7]);
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    colormap('parula');  
    title('Pupil Diameter');  
        set(gca, 'YTick', []); 

        
    subplot(1,5,4); phRasterFromTE(TE, all_behavior_trials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    set(gca, 'YTick', [], 'XLim', [-4 7]); 

    title('ChAT'); xlabel('Time frome odor (s)');
    subplot(1,5,5); phRasterFromTE(TE, all_behavior_trials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    set(gca, 'YTick', [], 'XLim', [-4 7]); 
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('DAT');    
    ax = findobj(gcf, 'Type', 'Axes');
    formatFigurePoster([12 4], [], 20);
        if saveOn    
        saveas(gcf, fullfile(savepath, savename), 'jpeg');
        saveas(gcf, fullfile(savepath, savename), 'epsc');
        saveas(gcf, fullfile(savepath, savename), 'fig');        
        end
%     set(ax, 'YLim', [40 80]);
%% allBehavior_averages_SFN
%     ylim = [-2 8];
    savename = 'allBehavior_averages';
    h=ensureFigure(saveName, 1); 
    
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-4 3], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording'};
    axh = [];
    subplot(1, 5, 1); [ha, hl] = plotEventAverageFromTE(TE, csPlusTrials & rewardTrials & expectTrials, 'Port1In', varargin{:});    

% [171 55 214]/256 [237 125 49]/256
    subplot(1, 5, 4);
    [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & expectTrials}, 1,...
    'FluorDataField', 'ZS', 'window', [-4, 3], 'cmap', [171 55 214]/256, 'alpha', 0); %high value, reward
     set(gca, 'XLim', [-4, 3]);%set(gca, 'YLim', ylim);


    subplot(1, 5, 5);
    [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & expectTrials}, 2,...
        'FluorDataField', 'ZS', 'window', [-4, 3], 'cmap', [237 125 49]/256, 'alpha', 0); %high value, reward
    legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    xlabel('time from cue (s)');      set(gca, 'XLim', [-4, 3]);

    formatFigurePoster([12 2], [], 20);
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end    
    %%
    ensureFigure('all_behavior_CsMinus', 1);
    reversals = find(diff(TE.BlockNumber(csMinusTrials, :))) + 1;
    subplot(1,4,1);
    title('Velocity');
    image(TE.Wheel.data.V(csMinusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
    set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 2, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 2]); 
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    subplot(1,4,2);
    image(TE.pupil.pupDiameterNorm(csMinusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
    set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    colormap('parula');  
    title('Pupil Diameter');    
    subplot(1,4,3); phRasterFromTE(TE, csMinusTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('ChAT'); xlabel('Time frome odor (s)');
    subplot(1,4,4); phRasterFromTE(TE, csMinusTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('DAT');    
    
    
    %% plot chat, dat, and pupil averages for cue responses
    
    savename = 'glm_averages';
    ensureFigure(savename, 1); ax = axes;
    yyaxis right;  ax.YColor = [0 0 0]; % ylabel('Pupil Diameter (norm.)');
    plotPupilAverageFromTE(TE, avgTrials & ~isnan(TE.pupil_cs), 'linespec', {'k'}, 'window', [-1 3]); hold on;
    set(gca, 'YLim', [1 1.12]);
    yyaxis left; ax.YColor = [0 0 0]; hold on;
    [ha, hl] = phPlotAverageFromTE(TE, glmTrials, 1,...
    'FluorDataField', 'ZS', 'window', [-1, 3], 'cmap', [171 55 214]/256, 'alpha', 0); %high value, reward
     set(gca, 'XLim', [-4, 3]);%set(gca, 'YLim', ylim);


    [ha, hl] = phPlotAverageFromTE(TE, glmTrials, 2,...
        'FluorDataField', 'ZS', 'window', [-1, 3], 'cmap', [237 125 49]/256, 'alpha', 0); %high value, reward
    xlabel('time from cue (s)');      set(gca, 'XLim', [-1, 3]);
    formatFigurePoster([5, 3], [], 24);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    
    %% make glm (NOT really a glm though, not using any linking functions)
    chat_bl = TE.phPeakMean_baseline(1).data(glmTrials);
    pup_bl = TE.pupilBaseline(glmTrials);
    wheel_bl = TE.wheelBaseline(glmTrials);
    pup_cs = TE.pupil_cs(glmTrials); % but mouse blinks sometimes so don't include this as regressor initially
    chat_cs = TE.phPeakMean_cs(1).data(glmTrials);
    dat_cs = TE.phPeakMean_cs(2).data(glmTrials);
    table_chat = table(pup_bl, wheel_bl, chat_cs);
    
    mdl_chat = fitglm(table_chat);
    
    %% ChAT Model
    glmTrials = csPlusTrials & rewardTrials & hitTrials;
    shiftPoints = 0.3 * 20;  % 0.3s lag and 20Hz sample rate
    pupScatter = TE.pupil.pupDiameter(:,1 + shiftPoints:end);
    chatScatter = TE.Photometry.data(1).ZS(:,1:end - shiftPoints);
    ensureFigure('pup_chat_scatter', 1); scatter(pupScatter(:),chatScatter(:), '.');
    
    glmEndTime = 3;

    input_pup = TE.pupil.pupDiameterNorm(glmTrials, bpX2pnt(-1, 20, -4) + shiftPoints :bpX2pnt(glmEndTime, 20, -4) + shiftPoints);
    input_pup = input_pup'; % just concatenate the trials together
    input_pup = input_pup(:); % just concatenate the trials together
    include = ~isnan(input_pup); 

    input = TE.Photometry.data(1).ZS(glmTrials, bpX2pnt(-1, 20, -4):bpX2pnt(glmEndTime, 20, -4));    
    input = input';
    numPoints = size(input, 1);
    numTrials = size(input, 2);
    x = diag(ones(1,numPoints));
    x2 = repmat(x, numTrials, 1);
    input = input(:);
     
    input = zscore(input(include));
    input_pup = zscore(input_pup(include));
    x2 = x2(include,:);

    % regression with design matrix (time as the regressor) 
    
    B_chat = regress(input,x2);
    ensureFigure('test', 1); plot(B_chat);
    
    
    % regression with time and pupil
    
    x3 = [x2 input_pup];
    
    B1_chat = regress(input,x3);
    ensureFigure('test', 1); plot(B1_chat);
    
    %
    shuff_pup = input_pup(randperm(length(input_pup)));
    shuff_x2 = x2(randperm(size(x2,1)),randperm(size(x2,2)));
    
    nrFolds = 10;
    for iFolds = 1:nrFolds
        idx = randperm(length(input));
        idx = idx(1:round(length(input)/10));
        cIdx = false(1,length(input));
        cIdx(idx) = true;
        
        B1_chat = regress(input(~cIdx),[x2(~cIdx,:) input_pup(~cIdx)]); %full model
        B2_chat = regress(input(~cIdx),[x2(~cIdx,:) shuff_pup(~cIdx)]); %time model
        B3_chat = regress(input(~cIdx),[shuff_x2(~cIdx,:) input_pup(~cIdx)]); %pupil model
        
        Y1{iFolds} = B1_chat' * ([x2(cIdx,:) input_pup(cIdx)])';
        Y2{iFolds} = B2_chat' * ([x2(cIdx,:) shuff_pup(cIdx)])';
        Y3{iFolds} = B3_chat' * ([shuff_x2(cIdx,:) input_pup(cIdx)])';
        Y{iFolds} = input(cIdx);
        
    end
    
    fullY = cat(1,Y{:});
    fullY1 = cat(2,Y1{:})';
    fullY2 = cat(2,Y2{:})';
    fullY3 = cat(2,Y3{:})';
    
    Rsq_chat(1) = corr2(fullY,fullY1).^2;
    Rsq_chat(2) = corr2(fullY,fullY2).^2;
    Rsq_chat(3) = corr2(fullY,fullY3).^2;
    
%% DAT Model
    shiftPoints = 0.3 * 20;  % 0.3s lag and 20Hz sample rate
    pupScatter = TE.pupil.pupDiameter(:,1 + shiftPoints:end);
    datScatter = TE.Photometry.data(2).ZS(:,1:end - shiftPoints);
    ensureFigure('pup_dat_scatter', 1); scatter(pupScatter(:),datScatter(:), '.');
    

    input = TE.Photometry.data(2).ZS(glmTrials, bpX2pnt(-1, 20, -4):bpX2pnt(glmEndTime, 20, -4));
    input_pup = TE.pupil.pupDiameterNorm(glmTrials, bpX2pnt(-1, 20, -4) + shiftPoints :bpX2pnt(glmEndTime, 20, -4) + shiftPoints);
    input_pup = input_pup';
    input_pup = input_pup(:);
    include = ~isnan(input_pup);
    
    input = input';
    numPoints = size(input, 1);
    numTrials = size(input, 2);
    x = diag(ones(1,numPoints));
    x2 = repmat(x, numTrials, 1);
    input = input(:);
     
    input = zscore(input(include));
    input_pup = zscore(input_pup(include));
    x2 = x2(include,:);

    % regression with design matrix (time as the regressor) 
    
    B_dat = regress(input,x2);
    ensureFigure('test', 1); plot(B_dat);
    
    
    % regression with time and pupil
    
    x3 = [x2 input_pup];
    
    B1_dat = regress(input,x3);
    ensureFigure('test', 1); plot(B1_dat);
    
    %
    shuff_pup = input_pup(randperm(length(input_pup)));
    shuff_x2 = x2(randperm(size(x2,1)),randperm(size(x2,2)));
    
    nrFolds = 10;
    for iFolds = 1:nrFolds
        idx = randperm(length(input));
        idx = idx(1:round(length(input)/10));
        cIdx = false(1,length(input));
        cIdx(idx) = true;
        
        B1_dat = regress(input(~cIdx),[x2(~cIdx,:) input_pup(~cIdx)]); %full model
        B2_dat = regress(input(~cIdx),[x2(~cIdx,:) shuff_pup(~cIdx)]); %time model
        B3_dat = regress(input(~cIdx),[shuff_x2(~cIdx,:) input_pup(~cIdx)]); %pupil model
        
        Y1{iFolds} = B1_dat' * ([x2(cIdx,:) input_pup(cIdx)])';
        Y2{iFolds} = B2_dat' * ([x2(cIdx,:) shuff_pup(cIdx)])';
        Y3{iFolds} = B3_dat' * ([shuff_x2(cIdx,:) input_pup(cIdx)])';
        Y{iFolds} = input(cIdx);
        
    end
    
    fullY = cat(1,Y{:});
    fullY1 = cat(2,Y1{:})';
    fullY2 = cat(2,Y2{:})';
    fullY3 = cat(2,Y3{:})';
    
    Rsq_dat(1) = corr2(fullY,fullY1).^2;
    Rsq_dat(2) = corr2(fullY,fullY2).^2;
    Rsq_dat(3) = corr2(fullY,fullY3).^2;    
    
    %% plot beta weights for ChAT and DAT
    
    
    xdata = [linspace(-1, 3, 81) 4];
    savename = 'glm_betaWeights_both';
    ensureFigure(savename, 1);
    subplot(1,2,1);
    hl = plot(xdata, B1_chat); hold on
    hl.Color = [171 55 214]/256; hl.LineWidth = 2;
        set(gca, 'XLim', [-1 4.5]);
%     formatFigurePoster([5, 3], [], 24);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end

%     savename = 'glm_betaWeights_dat';
%     ensureFigure(savename, 1);
subplot(1,2,2);
hl = plot(xdata, B1_dat); hold on
    hl.Color = [237 125 49]/256; hl.LineWidth = 2;
    set(gca, 'XLim', [-1 4.5]);
    formatFigurePoster([9, 3], [], 24);
    
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

%% bar graphs of R2 values

    savename = 'glm_Rsq_both';
    ensureFigure(savename, 1);
    subplot(1,2,1);
    bar(Rsq_chat, 'FaceColor', [171 55 214]/256);
    set(gca, 'YLim', [0 0.4], 'XTickLabel', {'full', 'time', 'pupil'});
    ylabel('R squared');
    subplot(1,2,2);
    bar(Rsq_dat, 'FaceColor', [237 125 49]/256);    
    set(gca, 'YLim', [0 0.4], 'XTickLabel', {'full', 'time', 'pupil'}, 'YTick', []);
    
    formatFigurePoster([9, 3], [], 24);
    
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    
   