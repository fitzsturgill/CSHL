% wheel_exampleMouse


DB = dbLoadExperiment('wheel');
saveOn = 1;


% load the session without reward

load('Z:\SummaryAnalyses\wheel_v1\DC_56_wheel_v1_Aug27_2018_Session1\TE.mat');

savepath = fullfile(DB.path, [filesep 'figure' filesep]);
ensureDirectory(savepath);



%% coherence
figSize = [2 1];

Photometry = 'PhotometryHF';


data_chat = TE.(Photometry).data(1).ZS';
% data_chat = diff(data_chat, 1, 1);
% data_chat = data_chat(:,2:end) - data_chat(:,1:end-1); % whiten
% data_chat = nanzscore2(data_chat); % standardize
% data_chat = nanzscore(data_chat); % standardize
data_dat = TE.(Photometry).data(2).ZS';
% data_dat = diff(data_dat, 1, 1);
% data_dat = data_dat(:,2:end) - data_dat(:,1:end-1); % whiten
% data_dat = nanzscore2(data_dat); % standardize
% data_dat = nanzscore(data_dat); % standardize


params.Fs = TE.(Photometry).sampleRate;
duration = length(TE.(Photometry).xData) / params.Fs;
params.trialave = 1;
params.err = [2 0.05];
params.tapers = [5 9];
params.pad = 1;
params.fpass = [1/duration 20];

[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_chat, data_dat, params);
f(1) = eps; % for log scale, you can't show zero
% scramble trial labels
si = randperm(size(data_dat, 2));

% [C_scram,phi_scram,S12_scram,S1_scram,S2_scram,f_scram,confC_scram, phistd_scram, Cerr_scram] = coherencyc(data_chat(:,si), data_dat, params);
% f_scram(1) = eps;

[C_scram,phi_scram,S12_scram,S1_scram,S2_scram,f_scram,confC_scram, phistd_scram, Cerr_scram] = coherencyc(circshift(data_chat, 1, 2), data_dat, params);
f_scram(1) = eps;

%
saveName = 'wheel_noRewards_coherence';
ensureFigure(saveName, 1);
axes; hold on;
boundedline(f_scram, C_scram, Cerr_scram(1,:)' - C_scram, 'k');%, 'alpha');
boundedline(f, C, Cerr(1,:)' - C, 'b');%, 'alpha');
% plot(f_scram, C_scram, 'k');%, 'alpha');
% plot(f_scram, Cerr_scram', 'Color', [0.8 0.8 0.8]);
% plot(f, C, 'b');%, 'alpha');
% plot(f, Cerr', 'Color', [0 0 0.8]);

set(gca, 'XScale', 'log', 'XLim', params.fpass);

xlabel('Frequency');
ylabel('Coherence');
%
formatFigurePublish('size', figSize);
if saveOn    
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));
end





















% 
% animal = 'DC_56';
% 
% 
% dbLoadAnimal(DB, animal);
% 
% 
% % cross correlogram from pre-reward period
% 
% 
% s1 = extractDataByTimeStamps(TE.Photometry.data(1).ZS, TE.Photometry.startTime, 20, TE.Reward, [-20 0]);
% s2 = extractDataByTimeStamps(TE.Photometry.data(2).ZS, TE.Photometry.startTime, 20, TE.Reward, [-20 0]);


