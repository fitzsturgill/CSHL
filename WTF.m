%% let's normalize us and cs responses to uncued reward, excluding ACh_17 for which there is basically no signal in one side of the brain
goodOnes = find(~ismember(DB.animals, 'ACh_17'));
nSub = length(goodOnes); 
nBoot = 1000;
alpha = 0.01;

% linecolors = [1 0 1; 0 1 0; 0 0 1; 1 0 0; 0 0.5 0; 0 0 0.5; 0 0 0.5];         
linecolors = [mycolors('puff'); mycolors('shock'); mycolors('reward'); mycolors('puff'); mycolors('shock'); mycolors('reward_cued'); mycolors('puff_cued')];
lineshapes = {'o'; 'o'; 's'; 's'; 's'; 'd'; 'd'};
trialSets = {'puff', 'shock', 'rew_cued', 'puff_cued', 'shock_cued', 'CSplus', 'CSminus_shock'};
trialSetNames = {'Air Puff', 'Shock', 'Reward cued', 'Puff cued', 'Shock cued', 'Cs+', 'CS-'};
% linecolors = [1 0 0; 0 0.5 0; 0 0 0.5; 0 0 0.5];         
% trialSets = {'puff_cued', 'shock_cued', 'CSplus', 'CSminus_shock'};
% trialSetNames = {'Puff cued', 'Shock cued', 'Cs+', 'CS-'};

allData = cell(nSub, length(trialSets));
allData_means = zeros(nSub, length(trialSets), 2);
allData_colors = repmat(linecolors, 1, 1, nSub);
allData_colors = permute(allData_colors, [3 1 2]); % put animal dimension first, then trial set dimension next, for sorting purposes later...
[allData_p, allData_h, allData_d] = deal(zeros(nSub, length(trialSets))); % p values, hypothesis, normalized auROC

allData_markers = repmat(lineshapes, 1, nSub);
allData_markers = permute(allData_markers, [2 1]);


denominator = us_pooled.rew.avg_mean(goodOnes,:);

% collect all the data
for acounter = 1:length(goodOnes)
    thisMouse = goodOnes(acounter);
    for counter = 1:length(trialSets)    
        trialSet = trialSets{counter};
        if ~ismember(trialSet, {'CSplus', 'CSminus_shock'})
            xData = us_pooled.(trialSets{counter}).mean{thisMouse}(:,1) ./ denominator(acounter,1);
            yData = us_pooled.(trialSets{counter}).mean{thisMouse}(:,2) ./ denominator(acounter,2);
            if ~ismember(trialSet, {'shock', 'shock_cued'})
                xData = xData .* 0;
                yData = yData .* 0;
            end
        else
            xData = cs_pooled.(trialSets{counter}).mean{thisMouse}(:,1) ./ denominator(acounter,1);
            yData = cs_pooled.(trialSets{counter}).mean{thisMouse}(:,2) ./ denominator(acounter,2);
            xData = xData .* 0;
            yData = yData .* 0;
        end
        allData_means(acounter, counter, 1) = nanmean(xData);
        allData_means(acounter, counter, 2) = nanmean(yData);
        allData{acounter, counter}(:,1) = xData;
        allData{acounter, counter}(:,2) = yData;
%         [D, P] = rocarea(xData, yData, 'boot', nBoot, 'scale');
        % altenatively, run a T test, and calculate Dprime yield effect
        % sizes
        D = abs(mean(xData) - mean(yData)) / sqrt((std(xData)^2 + std(yData)^2)/2);
        [~, P] = ttest(xData, yData);
        allData_p(acounter, counter) = P;
        allData_d(acounter, counter) = D;
    end
end

% let's do the Benjamini-Hochberg stepdown procedure to control false
% discovery rate to be < alpha

% rejected = slopes.auROC(slopes.p < 0.05);
% accepted = slopes.auROC(slopes.p >= 0.05);
[p_sorted, I] = sort(allData_p(:));
% auROC_sorted = auROC_all(I);
nComp = numel(allData_p);
alpha_adjusted = (1:nComp)' ./ nComp * alpha;
h_sorted = p_sorted >= alpha_adjusted;

allData_means = reshape(allData_means, [nComp, 2]);
allData_means = allData_means(I,:);
allData_colors = reshape(allData_colors, [nComp, 3]);
allData_colors = allData_colors(I,:);
allData_markers = reshape(allData_markers, [nComp, 1]);
allData_markers = allData_markers(I);
allData_d = allData_d(:);
allData_d = allData_d(I);