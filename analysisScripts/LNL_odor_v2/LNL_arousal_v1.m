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
    

    %
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

%% model with delay lines for pupil and odor delivery, god know if this is working correctly...
    modelTrials = csPlusTrials & rewardTrials & hitTrials;    
    tdlWindowOdor = [0 1]; % window for tapped delay line, window must be equal to or smaller than data window
    tdlWindowPupil = [-1 0];
    dataWindow = [-1 2]; % range of target data 0 1
    dataFullWindowOdor = tdlWindowOdor + dataWindow; % full window, allowing for delay lines
    dataFullWindowPupil = tdlWindowPupil + dataWindow; % full window, allowing for delay lines    

    % target
    input = TE.Photometry.data(1).ZS(modelTrials, bpX2pnt(dataWindow(1), 20, -4):bpX2pnt(dataWindow(2), 20, -4));    
    input = input';    
    numPoints = size(input, 1);
    numTrials = size(input, 2);   
    inputTemp = input;
    input = input(:);
%     input = zscore(input);
    
    % regressors
    input_pup = TE.pupil.pupDiameterNorm(modelTrials, bpX2pnt(dataFullWindowPupil(1), 20, -4):bpX2pnt(dataFullWindowPupil(2), 20, -4));
    input_pup = input_pup'; % nSamples x nTrials
%     input_pup = zeros(diff(dataFullWindowPupil) * 20, numTrials);
%     input_pup(bpX2pnt(0, 20, dataFullWindowPupil(1)),:) = 1; % odor time
    [pup_delayLine, truncPoints] = makeDelayLineDesignMatrix(input_pup, 'Fs', 20, 'window', tdlWindowPupil);
    
    input_odor = zeros(diff(dataFullWindowOdor) * 20 + 1, numTrials);
    input_odor(bpX2pnt(0, 20, dataFullWindowOdor(1)),:) = 1; % odor time
    [odor_delayLine, truncPoints] = makeDelayLineDesignMatrix(input_odor, 'Fs', 20, 'window', tdlWindowOdor);    
    
    yInt = ones(size(input));

    % regression with design matrix (time as the regressor) 
    ensureFigure('test', 1); subplot(1,2,1); imagesc(odor_delayLine(1:numPoints,:));
    subplot(1,2,2); imagesc(pup_delayLine(1:numPoints,:));
    
    %%
    X = [yInt odor_delayLine pup_delayLine];
    [B_chat bint, r, rint, stats] = regress(input, X);
    ensureFigure('test', 1); plot(B_chat);
    Rsq = regressCrossValidate(input, X);
    disp(['Rsq is ' num2str(Rsq)]);
    %%
    X = [yInt pup_delayLine];
    [B_chat bint, r, rint, stats] = regress(input, X);
    ensureFigure('test', 1); plot(B_chat, '-o');
    Rsq = regressCrossValidate(input, X);
    disp(['Rsq is ' num2str(Rsq)]);
    %%
    X = [yInt odor_delayLine];
    [B_chat bint, r, rint, stats] = regress(input, X);
    ensureFigure('test', 1); plot(B_chat);
    Rsq = regressCrossValidate(input, X);
    disp(['Rsq is ' num2str(Rsq)]);
    %%
    
    
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
    
    
    
    
    %% ChAT and DAT Models, fixed pupil lag
    modelTrials = csPlusTrials & rewardTrials & hitTrials;
    shiftPoints = 0.3 * 20;  % 0.3s lag and 20Hz sample rate
    
    modelEndTime = 3;

    input_pup = TE.pupil.pupDiameterNorm(modelTrials, bpX2pnt(-1, 20, -4) + shiftPoints :bpX2pnt(modelEndTime, 20, -4) + shiftPoints);
    input_pup = input_pup'; % just concatenate the trials together
    input_pup = input_pup(:); % just concatenate the trials together
    include = ~isnan(input_pup); 

    input = TE.Photometry.data(1).ZS(modelTrials, bpX2pnt(-1, 20, -4):bpX2pnt(modelEndTime, 20, -4));    
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
    
% DAT Model
    shiftPoints = 0.3 * 20;  % 0.3s lag and 20Hz sample rate
    pupScatter = TE.pupil.pupDiameter(:,1 + shiftPoints:end);
    datScatter = TE.Photometry.data(2).ZS(:,1:end - shiftPoints);
    ensureFigure('pup_dat_scatter', 1); scatter(pupScatter(:),datScatter(:), '.');
    

    input = TE.Photometry.data(2).ZS(modelTrials, bpX2pnt(-1, 20, -4):bpX2pnt(modelEndTime, 20, -4));
    input_pup = TE.pupil.pupDiameterNorm(modelTrials, bpX2pnt(-1, 20, -4) + shiftPoints :bpX2pnt(modelEndTime, 20, -4) + shiftPoints);
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
   