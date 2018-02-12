load('C:\Users\Adam\Dropbox\KepecsLab\_Fitz\PFF\AR.mat');

%% finding learning rate with correlation
% load('AR.mat');
load('C:\Users\Adam\Dropbox\KepecsLab\_Fitz\PFF\AR.mat');
param = KTD_defparam;
numCues = 2;
figure;
%SET DATA SOURCES
firstHalf = AR.allTrials;
secondHalf = AR.allTrials;
mainData = [firstHalf.csLicks.before secondHalf.csLicks.after];
reversalPoint = size(firstHalf.csLicks.before, 2);
numTrials = size(mainData,2);
%number of reversals
n = size(mainData,1);
rewards = [firstHalf.ReinforcementOutcome.before secondHalf.ReinforcementOutcome.after];
valves = [firstHalf.OdorValveIndex.before secondHalf.OdorValveIndex.after];
%SET VIEW RANGES FOR GRAPH
viewBefore = 10;
viewAfter = 20;
viewRange = (1:(viewAfter + viewBefore + 1))-(viewBefore + 1);
%ensureFigure('LearningRate',1);
colormap jet;
cmap = colormap;
%SETUP LEARNING RATES
lr = linspace(0.01,1,99);
corrs = zeros(length(lr),1);
for lrc = 1:length(lr)
    datas = nan(n,numTrials);
    models = nan(n,numTrials);
    zscores = zeros(n,numTrials);
    totalCorr = 0;
    for reversal = 1:n
        data = mainData(reversal,:);
        notnans = find(~isnan(data));
        data = data(notnans(1):notnans(end));
        numTrials = length(data); 
        rw = rewards(reversal,:);
        rw = rw(notnans(1):notnans(end));
        reward = zeros(numTrials,1);
        for i = (1:numTrials)
          if(rw(i) == "Reward")
               reward(i) = 1;
          elseif(rw(i) == "Punish")
                reward(i) = -1;
          else
             reward(i) = 0;
          end
        end
        X = zeros(numTrials,2);
        X(valves(reversal,notnans(1):notnans(end)) == 1,1) = 1;
        X(valves(reversal,notnans(1):notnans(end)) == 2,2) = 1;
        param.s = lr(1,lrc);
        param.q = 0.01;
        param.TD = 1;
        param.lr = lr(1,lrc);
        model = kalmanRW(X,reward,param);    
        rhat = zeros(n,1); % Predicted reward
        pe = zeros(numTrials,1); % prediction error
        w = zeros(numCues,numTrials); % weights
        Kn = zeros(numCues,numTrials); % Kalman gain
        offDiag = zeros(numTrials,1); % off diagonal term in posterior weight covariance matrix
        onDiag = zeros(2,numTrials); % on diagonal terms
        value = zeros(1,numTrials);
        for counter = 1:numTrials
           rhat(counter) = model(counter).rhat;
           pe(counter) = model(counter).dt;
           w(:,counter) = model(counter).w0;
           Kn(:,counter) = model(counter).K;
           offDiag(counter) = model(counter).C(2,1); % covariance matrix is symmetric so bottom left or top right corner of 2,2 matrix are equivalent
           onDiag(:,counter) = [model(counter).C(1,1); model(counter).C(2,2)];
           value(1,counter) = X(counter,:) * (model(counter).w0);   
        end
        models(reversal,notnans(1):notnans(end)) = value;
    end
    plotModel = nanmean(nanzscore((models(:,viewRange + reversalPoint)), 0, 2));
    plotData = nanmean(nanzscore((mainData(:,viewRange + reversalPoint)), 0, 2));
    
    plotModel = nanzscore((models(:,viewRange + reversalPoint)), 0, 2);
    plotData = nanzscore((mainData(:,viewRange + reversalPoint)), 0, 2);
    plotDiff = plotModel - plotData;
    absDiff = abs(plotDiff);
    sqrDiff = plotDiff .^2;
%     plotData = zscore(nanmean(mainData(:,viewRange + reversalPoint)));
%     corrs(lrc) = corr(plotModel', plotData');
    subplot(2,1,1);
    hold on;
    plot(viewRange, smoothdata(nanmean(sqrDiff),'movmedian',5),'Color',cmap((ceil((lrc/length(lr)) * 64)),:));
    hold off;
    subplot(2,1,2);
    hold on;
    plot(viewRange,plotModel, 'Color',cmap((ceil((lrc/length(lr)) * 64)),:),'LineWidth',.8);
    hold off;
    axis([-viewBefore viewAfter -2 2]);
end
%subplot(2,1,1);
%plot(lr, corrs);

subplot(2,1,2);
hold on;
plot(viewRange,plotData, 'Color','g','LineWidth',.8);
hold off;




%{

ensureFigure('Learning Rate Graph',1);
subplot(2,1,1);
hold on;
for i = 1:n
    xvals = (1:size(mainDataPlus,2)) - size(AR.csMinus.csLicks.before,2);
    plot(xvals(~isnan(mainDataPlus(i,:))),kalman_noise_filter(mainDataPlus(i,~isnan(mainDataPlus(i,:)))'), 'Color', 'g');
end
boundedline((1:size(mainDataPlus,2)) - size(AR.csMinus.csLicks.before,2) , nanmean(mainDataPlus), nansem(mainDataPlus), 'alpha');
ylabel("Anticipation (CS+)");
xlabel("Trials (Relative to Reversal)");
hold off;
subplot(2,1,2);
hold on;
for i = 1:size(mainDataMinus,1)
    xvals = (1:size(mainDataMinus,2)) - size(AR.csPlus.csLicks.before,2);
    plot(xvals(~isnan(mainDataMinus(i,:))),kalman_noise_filter(mainDataMinus(i,~isnan(mainDataMinus(i,:)))'), 'Color', 'g');
end
boundedline((1:size(mainDataMinus,2)) - size(AR.csPlus.csLicks.before,2) , nanmean(mainDataMinus), nansem(mainDataMinus), 'alpha'); 
ylabel("Anticipation (CS-)");
xlabel("Trials (Relative to Reversal)");
hold off;


lastCorr;
ensureFigure('learning Rate');
plot(corrs);
%}