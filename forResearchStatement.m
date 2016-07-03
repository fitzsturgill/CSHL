
%     clear

%     cd('/Users/fitz/Documents/KepecsLab/BehaviorData/Chat_GCaMP_pink/SO_Training_NIDAQ/Session Data');
%     pinkSessions(1) = load('Chat_GCaMP_pink_SO_Training_NIDAQ_Sep25_2015_Session2.mat');
%     pinkSessions(2) = load('Chat_GCaMP_pink_SO_Training_NIDAQ_Sep28_2015_Session2.mat');  
%     pinkSessions(3) = load('Chat_GCaMP_pink_SO_Training_NIDAQ_Sep29_2015_Session1.mat');

%% First GCaMP
%     sessionTitle = 'pinkSessions';
%     fig = figure('Name', sessionTitle);
%     title_ax = textAxes(fig, sessionTitle, 12); %fontsize = 12
%     
%     stimBarY = -0.005;
%     yLims = [-0.01 0.01];
% 
%     %primary comparison: toneA reward vs toneB punish
% 
%     ax(1) = axes; hold on
%     [avg_g avgX_g avgSEM_g allTrials_g ] = sessionsAvg(pinkSessions, 1, 2, 'g', ax(1));
%     [avg_r avgX_r avgSEM_r allTrials_r ] = sessionsAvg(pinkSessions, 2, 4, 'r', ax(1));
% %     textBox('green: toneA, reward; red: toneB, punish', [],[], 12);
%     xlabel('time (s)');
%     ylabel('GCaMP norm');
%     addStimulusBar(ax(1), [0 1 stimBarY], '', [0 0 0], 3);
%     addStimulusBar(ax(1), [1.5 1.7 stimBarY], '', [0 0 0], 3);
    
    %% Second: Licks

    % I hard coded 0.1sec as my sample size for the lick histogram counts
    [rewardLickTimes rewardLickTrials rewardLickRate_avg rewardCountsX rewardLickRate_sem] = researchStatementLicks(pinkSessions, 1, 2);
    [punishLickTimes punishLickTrials punishLickRate_avg punishCountsX punishLickRate_sem] = researchStatementLicks(pinkSessions, 2, 4);
    
    
    

    figure; 
    ax=axes(...
    'YDir', 'reverse'...
    );
    linecustommarker(rewardLickTimes, rewardLickTrials, [], [], ax);  
    title('reward');

    figure; 
    ax=axes(...
    'YDir', 'reverse'...
    );
    linecustommarker(punishLickTimes, punishLickTrials, [], [], ax);  
    title('punish');
% 
%     pCounts = histcounts(punishLickTimes, -2:0.5:7);
%     rCounts = histcounts(rewardLickTimes, -2:0.5:7);
%     figure; plot(-2:0.5:6.9, pCounts); hold on
%     plot(-2:0.1:6.5, rCounts); 
    
    
    % make trial histogram
    figure; 
    ax = axes; hold on
    plot(rewardCountsX, rewardLickRate_avg, 'g');
    plot(punishCountsX, punishLickRate_avg, 'r');
%     boundedline(rewardCountsX, rewardLickRate_avg, rewardLickRate_sem, 'g', ax, 'alpha');
%     boundedline(punishCountsX, punishLickRate_avg, punishLickRate_sem, 'r', ax, 'alpha');
%     figure; plot(rewardCountsX, mean(rewardCounts));
    
%     figure; plot(punishCountsX, mean(punishCounts));




        
%     linecustommarker(lickTimes, lickTrials, [], [], ax);  
%     title(label);
        
