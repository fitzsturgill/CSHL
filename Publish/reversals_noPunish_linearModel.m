DB = dbLoadExperiment('reversals_noPunish_publish');
saveOn = 1;
savePath = fullfile(DB.path, 'pooled', filesep);
ensureDirectory(savePath);
animals = {'DC_44'; 'DC_46'; 'DC_47'; 'DC_53'; 'DC_54'; 'DC_56'}; 

nm = length(animals);
%
window = [-1 2];
Fs = 20;
kernalSize = diff(window) * Fs + 1;
rsq_all = zeros(nm, 3, 2);
[bFull_all, bPupil_all, bTime_all] = deal(zeros(nm, kernalSize*2, 2));
% [bPupil_all bTime_all] = deal(zeros(nm, kernalSize, 2));


winIx = bpX2pnt(window(1), 20, -4):bpX2pnt(window(2), 20, -4);
for counter = 1:nm
    animal = animals{counter};
    dbLoadAnimal(DB, animal);
    if strcmp(animal, 'DC_53')
        lmTrials = csPlusTrials & rewardTrials & hitTrials & (TE.sessionIndex ~= 4); % exclude session 4 where mouse runs for entire session, k
    else
        lmTrials = csPlusTrials & rewardTrials & hitTrials;
    end
%     figure(hp);
%     subplot(2,3,counter); hold on;
%     plotPupilAverageFromTE(TE, lmTrials & ~any(isnan(TE.pupil.pupDiameterNorm), 2), 'linespec', {'k'}, 'window', [-3.9 6.9]);
    
    


    m = sum(lmTrials);
    n = length(winIx);
%     shiftPoints = 0.3 * 20;  % 0.3s lag and 20Hz sample rate
    for ch = 1:2



        % regressors
        
        % pupil 
        shift = floor(length(winIx)/2); 
        shiftIx = -shift:shift;
        o = length(shiftIx);
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
        pup_shuff = pup(randperm(numTrials), randperm(numPoints));
        
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

%%
hf = ensureFigure('kernel_full');
hp = ensureFigure('kernel_pupil');
ht = ensureFigure('kernel_time');
hr = ensureFigure('Rsq_byMouse');
ch = 1;
for counter = 1:nm
    figure(hf);
    subplot(2,3,counter); hold on;
    plot(bFull_all(counter, 1:n, ch));
    plot(bFull_all(counter, n+1:end, ch));
    textBox(animals{counter});
    
    figure(hp);
    subplot(2,3,counter); hold on;
    plot(bPupil_all(counter, 1:n, ch), '.');
    plot(bPupil_all(counter, n+1:end, ch));
    textBox(animals{counter});
    
    figure(ht);
    subplot(2,3,counter); hold on;
    plot(bTime_all(counter, 1:n, ch));
    plot(bTime_all(counter, n+1:end, ch), '.');    
    textBox(animals{counter});
    
    figure(hr);
    subplot(2,3,counter); hold on;  
    bar(squeeze(rsq_all(counter, :,1)));
    set(gca, 'XTickLabel', {'full', 'time', 'pup'}, 'XTick', [1 2 3]);
    textBox(animals{counter});
end

%% bar graphs
ensureFigure('Rsquared_ch1');
plot(squeeze(rsq_all(:,:,1))', 'o-');

ensureFigure('Rsquared_ch2');
plot(squeeze(rsq_all(:,:,2))', 'o-');
    
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