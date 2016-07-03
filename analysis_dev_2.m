% load script


%     clear
%% switch directory, gather experiment names
    cd('/Users/fitz/Documents/KepecsLab/BehaviorData/Chat_GCaMP_pink/SO_Training_NIDAQ/Session Data');
    clear sessionNames
    sessionNames{1} = 'Chat_GCaMP_pink_SO_Training_NIDAQ_Sep25_2015_Session2.mat';
    sessionNames{2} = 'Chat_GCaMP_pink_SO_Training_NIDAQ_Sep28_2015_Session2.mat';
    sessionNames{3} = 'Chat_GCaMP_pink_SO_Training_NIDAQ_Sep29_2015_Session1.mat';
    sessionNames{4} = 'Chat_GCaMP_pink_SO_Training_NIDAQ_Oct02_2015_Session1.mat';
    
%%  
    session = 1;
    sessionName = sessionNames{session};
    load(sessionName);
    pavlovian_lick_rasters(SessionData, sessionName);
    pavlovian_Photometry(SessionData, sessionName);