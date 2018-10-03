sessions = bpLoadSessions([], 'Dummy Subject_wheel_v1_Sep24_2018_Session4.mat', 'Z:\FitzRig2\Data\Dummy Subject\wheel_v1\Session Data');

%% 
TE = makeTE_wheel_v1(sessions); % make TE structure

channels=[]; dFFMode = {}; %BL = {};
% if sessions(1).SessionData.Settings.GUI.LED1_amp > 0
    channels(end+1) = 1;
    dFFMode{end+1} = 'simple';
%     BL{end + 1} = [1 4];
% end

% if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple';
%     BL{end + 1} = [2 4];    
% end

baselineEnd = sessions(1).SessionData.NidaqData{1,2}.duration(1) - 0.2;
        
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', channels, 'baseline', [1 baselineEnd], 'startField', 'Baseline', 'downsample', 10, 'ACfilter', 1);%, 'forceAmp', 1);

TE.PhotometryDC = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', channels, 'baseline', [1 baselineEnd], 'startField', 'Baseline', 'downsample', 10, 'ACfilter', 0);
Fs_raw = 6100;


%% plot curve for Photometry DC
ensureFigure('DC_noise');
trials = [22 2:2:20];
noise_dc = zeros(length(trials), 2);
baseline_dc = zeros(length(trials), 2);
for counter = 1:length(trials)
    trial = trials(counter);
    noise_dc(counter, 1) = std(TE.PhotometryDC.data(1).raw(trial, :));
    noise_dc(counter, 2) = std(TE.PhotometryDC.data(2).raw(trial, :));
    baseline_dc(counter, 1) = mean(TE.PhotometryDC.data(1).raw(trial, :));
    baseline_dc(counter, 2) = mean(TE.PhotometryDC.data(2).raw(trial, :));
end
baseline_dc2 = baseline_dc(2:end, :) - baseline_dc(1,:); % subtract off dark volrate
noise_dc2 = noise_dc(2:end,:); % doesn't include LED off
plot(baseline_dc, noise_dc, '+-');

% ensureFigure('justthetrials', 1); plot(TE.PhotometryDC.data(2).raw(trials, :)');
    
    
%% plot curve for Photometry AC
ensureFigure('AC_noise');
trials = [22 2:2:20] - 1;
noise_ac = zeros(length(trials), 2);
baseline_ac = zeros(length(trials), 2);
for counter = 1:length(trials)
    trial = trials(counter);
    noise_ac(counter, 1) = std(TE.Photometry.data(1).raw(trial, 610:end-610));
    noise_ac(counter, 2) = std(TE.Photometry.data(2).raw(trial, 610:end-610));
    baseline_ac(counter, 1) = mean(TE.Photometry.data(1).raw(trial, 610:end-610));
    baseline_ac(counter, 2) = mean(TE.Photometry.data(2).raw(trial, 610:end-610));
end
plot(baseline_ac, noise_ac);

ensureFigure('justthetrials', 1); plot(TE.Photometry.data(2).raw(trials, :)');

%% plot both AC and DC
 ensureFigure('both', 1); hold on;
 subplot(1,2,1);
 plot(baseline_dc2(:,1), noise_dc2(:,1), 'g'); hold on;
 plot(baseline_ac(2:end,1), noise_ac(2:end,1), 'k');
 legend('dc', 'ac'); xlabel('baseline (mean)'); ylabel('noise (std)');
 subplot(1,2,2);
 plot(baseline_dc2(:,2), noise_dc2(:,2), 'r'); hold on;
 plot(baseline_ac(2:end,2), noise_ac(2:end,2), 'k');
 legend('dc', 'ac'); xlabel('baseline (mean)'); ylabel('noise (std)');
 
 %% compare DC and AC single trials
 ensureFigure('compare', 1); plot(TE.Photometry.xData, TE.Photometry.data(1).raw(19, :), 'b'); hold on;
 plot(TE.Photometry.xData, TE.PhotometryDC.data(1).raw(20, :) - mean(TE.PhotometryDC.data(1).raw(22, :)), 'g'); 
 
 %% plot histograms, mean centered
 ensureFigure('compare2', 1); hold on;
 histogram(TE.Photometry.data(1).raw(19, 610:end-610) - mean(TE.Photometry.data(1).raw(19, 610:end-610)), linspace(-0.01, .01, 200));
 histogram(TE.PhotometryDC.data(1).raw(20, 610:end-610) - mean(TE.PhotometryDC.data(1).raw(20, 610:end-610)), linspace(-0.01, .01, 200)); 
 set(gca', 'XLim', [-0.005 0.005]);
 legend(sprintf('ACsigma=%.3g', std(TE.Photometry.data(1).raw(19, :))), sprintf('DCsigma=%.3g', std(TE.PhotometryDC.data(1).raw(20, :))));