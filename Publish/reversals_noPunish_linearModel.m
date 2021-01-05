DB = dbLoadExperiment('reversals_noPunish_publish');
saveOn = 1;
savePath = fullfile(DB.path, 'pooled', filesep);
ensureDirectory(savePath);
animals = {'DC_44'; 'DC_46'; 'DC_47'; 'DC_53'; 'DC_54'; 'DC_56'}; 

nm = length(animals);
%
window = [-1 2];

Fs = 20;

%% make histograms of hit trial lick counts for csPlus Odor

saveName = 'hitLick_histogram';
ensureFigure(saveName);
ha = [];
for counter = 1:nm
    
    animal = animals{counter};
    dbLoadAnimal(DB, animal);
    ha(end+1) = subplot(2,3,counter); hold on; 
    lickData = TE.licks_cs.count(hitTrials & csPlusTrials);
    tercile = percentile(lickData, 1/3);
    histogram(lickData);
    ylim = get(gca, 'YLim');
    line([tercile tercile], ylim, 'Color', 'k', 'LineWidth', 1.5)
    textBox(animal);
end
sameXScale(ha);

%%
kernalSize = diff(window) * Fs + 1;
rsq_all = zeros(nm, 3, 2);
binSize = 1;
winIx = bpX2pnt(window(1), 20, -4):bpX2pnt(window(2), 20, -4);
shift = floor(length(winIx)/2); 
shiftIx = -shift:binSize:shift;
o = length(shiftIx);
[bFull_all, bPupil_all, bTime_all] = deal(zeros(nm, kernalSize + o, 2));
% [bPupil_all bTime_all] = deal(zeros(nm, kernalSize, 2));
n = length(winIx);

%%
for counter = 1:nm
    animal = animals{counter};
    dbLoadAnimal(DB, animal);
    
    lickData = TE.licks_cs.count(csPlusTrials & hitTrials);
    tercile = percentile(lickData, 1/3);
    if strcmp(animal, 'DC_53')
        lmTrials = csPlusTrials & hitTrials & (TE.sessionIndex ~= 4) & (TE.licks_cs.count > tercile); % exclude session 4 where mouse runs for entire session, k
    else
        lmTrials = csPlusTrials & hitTrials & (TE.licks_cs.count > tercile);
    end
%     figure(hp);
%     subplot(2,3,counter); hold on;
%     plotPupilAverageFromTE(TE, lmTrials & ~any(isnan(TE.pupil.pupDiameterNorm), 2), 'linespec', {'k'}, 'window', [-3.9 6.9]);
    
    


    m = sum(lmTrials);

%     shiftPoints = 0.3 * 20;  % 0.3s lag and 20Hz sample rate
    for ch = 1:2



        % regressors
        
        % pupil 


        winIx2 = repmat(winIx, length(shiftIx), 1) + shiftIx(:);
        winIx2 = winIx2'; % put each copy of trial along columns, each copy is shifted by an amount mapped onto column number
        winIx2 = winIx2(:); % concatenate them
        pup = TE.pupil.pupDiameterNorm(lmTrials, winIx2);
        pup = reshape(pup, [m n o]);
        pup = permute(pup, [2 1 3]);
        pup = reshape(pup, [m*n, o]); % now you should have versions of concatenated pupil data shifted by different amounts along columns
        include = ~any(isnan(pup), 2);
        pup = pup(include,:);
        
        
        % response variable
        input = TE.PhotometryExpFit.data(ch).ZS(lmTrials, winIx);
        input = input';
        input = input(:); % concatenate trials together        
        input = input(include);
        numPoints = n;
        numTrials = sum(include);        

        
        % time
        x = diag(ones(1,numPoints));
        x2 = repmat(x, m, 1);
        x2 = x2(include, :);


        % zscore
        input = zscore(input);
        pup = zscore(pup);
        
        % make shuffled versions
        x2_shuff = x2(randperm(numTrials), randperm(numPoints));
        pup_shuff = pup(randperm(numTrials), randperm(size(pup, 2)));
        
        % regression full
        B_chat = regress(input,[x2 pup]);

        % regression with time and pupil
        B_time = regress(input,[x2 pup_shuff]);
        
        % regression with pupil
        B_pupil = regress(input,[x2_shuff pup]);        

%         figure(hk);
%         subplot(2,3,counter); hold on;
%         plot(B1_chat(1:n));  plot(B1_chat(n+1:end), 'o');
%         textBox(animal);

        % cross validation
        nrFolds = 10;
        for iFolds = 1:nrFolds
            idx = randperm(length(input));
            idx = idx(1:round(length(input)/nrFolds));
            cIdx = false(1,length(input));
            cIdx(idx) = true; % indexes for 1/nrFolds fraction of random time points

            % regression model on most of the data (nrFolds-1 / nrFolds)
            B1_chat = regress(input(~cIdx),[x2(~cIdx,:) pup(~cIdx, :)]); %full model
            B2_chat = regress(input(~cIdx),[x2(~cIdx,:) pup_shuff(~cIdx, :)]); %time model
            B3_chat = regress(input(~cIdx),[x2_shuff(~cIdx,:) pup(~cIdx, :)]); %pupil model

            Y1{iFolds} = B1_chat' * ([x2(cIdx,:) pup(cIdx, :)])';
            Y2{iFolds} = B2_chat' * ([x2(cIdx,:) pup_shuff(cIdx, :)])';
            Y3{iFolds} = B3_chat' * ([x2_shuff(cIdx,:) pup(cIdx, :)])';
            Y{iFolds} = input(cIdx);

        end

        fullY = cat(1,Y{:});
        fullY1 = cat(2,Y1{:})';
        fullY2 = cat(2,Y2{:})';
        fullY3 = cat(2,Y3{:})';

        Rsq(1) = corr2(fullY,fullY1).^2;
        Rsq(2) = corr2(fullY,fullY2).^2;
        Rsq(3) = corr2(fullY,fullY3).^2;
        rsq_all(counter, :, ch) = Rsq;
        bFull_all(counter, :, ch) = B_chat;
        bPupil_all(counter, :, ch) = B_pupil;
        bTime_all(counter, :, ch) = B_time;
    end
end

save(fullfile(savePath, 'linearModel.mat'), 'rsq_all', 'bFull_all', 'bPupil_all', 'bTime_all'); 



%% load the data, if desired
load(fullfile(savePath, 'linearModel.mat')); 

%%
hf = ensureFigure('kernel_full');
hp = ensureFigure('kernel_pupil');
ht = ensureFigure('kernel_time');
hr = ensureFigure('Rsq_byMouse');

for counter = 1:nm
    for ch = 1:2
        figure(hf);
        subplot(2,6,counter + 6*(ch-1)); hold on;
        plot(bFull_all(counter, 1:n, ch));
        plot(bFull_all(counter, n+1:end, ch));
        textBox(animals{counter});
        if counter == 1
            ylabel(['Channel ' num2str(ch)]);
        end
        
        figure(hp);
        subplot(2,6,counter + 6*(ch-1)); hold on;
        plot(bPupil_all(counter, 1:n, ch), '.');
        plot(bPupil_all(counter, n+1:end, ch));
        textBox(animals{counter});
        if counter == 1
            ylabel(['Channel ' num2str(ch)]);
        end
        
        figure(ht);
        subplot(2,6,counter + 6*(ch-1)); hold on;
        plot(bTime_all(counter, 1:n, ch));
        plot(bTime_all(counter, n+1:end, ch), '.');    
        textBox(animals{counter});
        if counter == 1
            ylabel(['Channel ' num2str(ch)]);
        end
        
        figure(hr);
        subplot(2,6,counter + 6*(ch-1)); hold on;  
        bar(squeeze(rsq_all(counter, :,ch)));
        set(gca, 'XTickLabel', {'full', 'time', 'pup'}, 'XTick', [1 2 3]);
        textBox(animals{counter});
        if counter == 1
            ylabel(['Channel ' num2str(ch)]);
        end        
    end
end
sameYScale(get(hf, 'Children'));
sameYScale(get(hp, 'Children'));
sameYScale(get(ht, 'Children'));
sameYScale(get(hr, 'Children'));

%% bar graphs
ensureFigure('Rsquared_ch1');
plot(squeeze(rsq_all(:,:,1))', 'o-');

ensureFigure('Rsquared_ch2');
plot(squeeze(rsq_all(:,:,2))', 'o-');

%% summary bar graphs

saveName = 'Rsquared_all';
ensureFigure(saveName);
rsq_SEM = squeeze(std(rsq_all) / sqrt(size(rsq_all, 1)));
rsq_mean = squeeze(mean(rsq_all));
% rsq_X = repmat((1:3)', 1, 2) + [-0.25 0.25];
% bar(rsq_X, rsq_mean); hold on;
% errorbar(rsq_X, rsq_mean, rsq_SEM);

h(1) = subplot(1,2,1);
barwitherr(rsq_SEM(:,1), rsq_mean(:,1), 'FaceColor', mycolors('chat')); hold on;
set(gca,'XTickLabel',{'Full','Time','Pupil'})
ylabel('Variance Explained');
textBox('ACh');
h(2) = subplot(1,2,2);
barwitherr(rsq_SEM(:,2), rsq_mean(:,2), 'FaceColor', mycolors('dat'));
set(gca,'XTickLabel',{'Full','Time','Pupil'})
textBox('DA');
sameYScale(h);

formatFigurePublish('size', [4 2]);
if saveOn    
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));
end

%% ChAT
% statistics, 1 factor ANOVA follwed by multiple comparisons, treat DA and
% ChAT seperately because hard to interpret interation, meaning that means
% of DA aren't necessarily a scaled version of ACh means

[anova1_chat,~, stats_chat] = anova1(squeeze(rsq_all(:,:,1)));
c = multcompare(stats_chat); % tukey-kramer, honest significant difference
colNames = {'group1', 'group2', 'Conf_lower', 'estimate', 'Conf_upper', 'p_value'};
c_chat = table(c(:,1),c(:,2),c(:,3),c(:,4),c(:,5),c(:,6), 'VariableNames', colNames);
%c is Matrix of multiple comparison results, returned as an p-by-6 matrix of 
% scalar values, where p is the number of pairs of groups. Each row of the 
% matrix contains the result of one paired comparison test. 
% Columns 1 and 2 contain the indices of the two samples being compared. 
% Column 3 contains the lower confidence interval, column 4 contains the estimate, 
% and column 5 contains the upper confidence interval. Column 6 contains the 
% p-value for the hypothesis test that the corresponding mean difference is not equal to 0.


%% DAT
[anova1_dat,~, stats_dat] = anova1(squeeze(rsq_all(:,:,2)));
c = multcompare(stats_dat); % tukey-kramer, honest significant difference
colNames = {'group1', 'group2', 'Conf_lower', 'estimate', 'Conf_upper', 'p_value'};
c_dat = table(c(:,1),c(:,2),c(:,3),c(:,4),c(:,5),c(:,6), 'VariableNames', colNames);





% statistics, 2 fator ANOVA followed by multiple comparisons,
% rsq_2factor = reshape(permute(rsq_all, [1 3 2]), 12, 3); % columns are different models, rows are ChAT vs Dop., 6 replicates
% [p, tbl, stats] = anova2(rsq_2factor, 6);
% follow up with multiple comparisons, rows and columns
% [c,m,h,nms] = multcompare(stats);
%% example mouse beta coefficients

        saveName = 'beta_Coefficients_example';
        ensureFigure(saveName);
        ha=zeros(6,1);
        mouse = 'DC_44';
        counter = find(strcmp(mouse, animals)); % mouse
        ha(1) = subplot(2,3,1);hold on;
        plot(bFull_all(counter, 1:n, 1), 'Color', mycolors('chat')); hold on;
        plot(bFull_all(counter, n+1:end, 1), 'Color', mycolors('chat'), 'LineWidth', 2);        
        ylabel('Beta coefficient');
        title('Full');
        
        ha(2) = subplot(2,3,2); hold on;
        plot(bTime_all(counter, 1:n, 1), 'Color', mycolors('chat'));
        plot(bTime_all(counter, n+1:end, 1), 'Color', mycolors('chat'), 'LineWidth', 2);        
        title('Time only (pup shuff)');
        
        ha(3) = subplot(2,3,3);hold on;
        plot(bPupil_all(counter, 1:n, 1), 'Color', mycolors('chat'));
        plot(bPupil_all(counter, n+1:end, 1), 'Color', mycolors('chat'), 'LineWidth', 2);                
        title('Pupil only (time shuff)');
        
        ha(4) = subplot(2,3,4);hold on;
        plot(bFull_all(counter, 1:n, 2), 'Color', mycolors('dat'));
        plot(bFull_all(counter, n+1:end, 2), 'Color', mycolors('dat'), 'LineWidth', 2);        
        ylabel('Beta coefficient');
        
        ha(5) = subplot(2,3,5);hold on;
        plot(bTime_all(counter, 1:n, 2), 'Color', mycolors('dat'));
        plot(bTime_all(counter, n+1:end, 2), 'Color', mycolors('dat'), 'LineWidth', 2);        
        
        ha(6) = subplot(2,3,6);hold on;
        plot(bPupil_all(counter, 1:n, 2), 'Color', mycolors('dat'));
        plot(bPupil_all(counter, n+1:end, 2), 'Color', mycolors('dat'), 'LineWidth', 2);                
        sameYScale(ha);
        formatFigurePublish('size', [8 4]);

        if saveOn    
            print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
            saveas(gcf, fullfile(savePath, [saveName '.fig']));
            saveas(gcf, fullfile(savePath, [saveName '.jpg']));
        end
        
        
   
%% example coefficients for figure
        saveName = 'beta_Coefficients_Figure';
        ensureFigure(saveName);
        hat=[];
        hap=[];
        xpp = size(bFull_all, 2) - n;
        xd1 = (0:(n-1)) ./ Fs;
        xd2 = diff(window)/2 - ((0:(xpp-1))./Fs);        
        xd1 = xd1 - 1;
        mouse = 'DC_44';
        counter = find(strcmp(mouse, animals)); % mouse
        hat(end+1) = subplot(2,2,1);hold on;
        plot(xd1, bFull_all(counter, 1:n, 1), 'Color', mycolors('chat')); hold on;
        set(gca, 'XTickLabel', {}, 'XTick', [-1 0 1 2]);
%         ylabel('ACh weights');
%         title('Time');
        
        hap(end+1) = subplot(2,2,2); hold on;
        plot(xd2, bFull_all(counter, n+1:end, 1), 'Color', mycolors('chat'));        
        set(gca, 'XTickLabel', {}, 'XTick', [-1 0 1]);
%         title('Pupil');
        
        hat(end+1) = subplot(2,2,3);hold on;
        plot(xd1, bFull_all(counter, 1:n, 2), 'Color', mycolors('dat')); hold on;
        set(gca,'XTick', [-1 0 1 2]);
%         ylabel('DA weights');
        
        hap(end+1) = subplot(2,2,4); hold on;
        plot(xd2, bFull_all(counter, n+1:end, 2), 'Color', mycolors('dat'));                
        set(gca, 'XTick', [-1 0 1]);
%         sameYScale(hat);sameYScale(hap);
        sameYScale([hat; hap]);
        formatFigurePublish('size', [3.35 2]);

        if saveOn    
            print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
            saveas(gcf, fullfile(savePath, [saveName '.fig']));
            saveas(gcf, fullfile(savePath, [saveName '.jpg']));
        end

    
%% Old version where I have a fully shifted version of pupil, very correlated regressors
%     %% ChAT Model
%     window = [-1 1];
%     winIx = bpX2pnt(window(1), 20, -4):bpX2pnt(window(2), 20, -4);
%     m = sum(lmTrials);
%     n = length(winIx);
% %     shiftPoints = 0.3 * 20;  % 0.3s lag and 20Hz sample rate
%     input = TE.Photometry.data(1).ZS(lmTrials, winIx);
% %     input_pup = TE.pupil.pupDiameterNorm(lmTrials, winIx + shiftPoints);
%     shift = floor(length(winIx)/2);
% 
%     shiftIx = -shift:20:shift;
%     o = length(shiftIx);
%     winIx2 = repmat(winIx, length(shiftIx), 1) + shiftIx(:);
%     winIx2 = winIx2'; % put each copy of trial along columns, each copy is shifted by an amount mapped onto column number
%     winIx2 = winIx2(:); % concatenate them
%     input_pup = TE.pupil.pupDiameterNorm(lmTrials, winIx2);
%     input_pup = reshape(input_pup, [m n o]);
%     input_pup = permute(input_pup, [2 1 3]);
%     input_pup = reshape(input_pup, [m*n, o]); % now you should have versions of concatenated pupil data shifted by different amounts along columns
%     
%     
%     
% %     input_pup = input_pup';
% %     input_pup = input_pup(:);
%     include = ~any(isnan(input_pup), 2);
%     
%     input = input';
%     numPoints = size(input, 1);
%     numTrials = size(input, 2);
%     x = diag(ones(1,numPoints));
%     x2 = repmat(x, numTrials, 1);
%     input = input(:);
% %   
%     input = zscore(input(include));
%     input_pup = zscore(input_pup(include, :));
%     x2 = x2(include,:);
% 
%     % regression with design matrix (time as the regressor) 
%     
%     B_chat = regress(input,x2);
%     ensureFigure('test', 1); plot(B_chat);
%     
%     
%     % regression with time and pupil
%     
%     x3 = [x2 input_pup];
%     
%     B1_chat = regress(input,x3);
%     ensureFigure('test', 1); plot(B1_chat(1:n)); hold on; plot(B1_chat(n+1:end));
%     
%     %
%     shuff_pup = input_pup(randperm(size(input_pup, 1)), :);
%     shuff_x2 = x2(randperm(size(x2,1)),randperm(size(x2,2)));
%     
%     nrFolds = 10;
%     for iFolds = 1:nrFolds
%         idx = randperm(length(input));
%         idx = idx(1:round(length(input)/nrFolds));
%         cIdx = false(1,length(input));
%         cIdx(idx) = true; % indexes for 1/nrFolds fraction of random time points
%         
%         % regression model on most of the data (nrFolds-1 / nrFolds)
%         B1_chat = regress(input(~cIdx),[x2(~cIdx,:) input_pup(~cIdx, :)]); %full model
%         B2_chat = regress(input(~cIdx),[x2(~cIdx,:) shuff_pup(~cIdx, :)]); %time model
%         B3_chat = regress(input(~cIdx),[shuff_x2(~cIdx,:) input_pup(~cIdx, :)]); %pupil model
%         
%         Y1{iFolds} = B1_chat' * ([x2(cIdx,:) input_pup(cIdx, :)])';
%         Y2{iFolds} = B2_chat' * ([x2(cIdx,:) shuff_pup(cIdx, :)])';
%         Y3{iFolds} = B3_chat' * ([shuff_x2(cIdx,:) input_pup(cIdx, :)])';
%         Y{iFolds} = input(cIdx);
%         
%     end
%     
%     fullY = cat(1,Y{:});
%     fullY1 = cat(2,Y1{:})';
%     fullY2 = cat(2,Y2{:})';
%     fullY3 = cat(2,Y3{:})';
%     
%     Rsq(1) = corr2(fullY,fullY1).^2;
%     Rsq(2) = corr2(fullY,fullY2).^2;
%     Rsq(3) = corr2(fullY,fullY3).^2;