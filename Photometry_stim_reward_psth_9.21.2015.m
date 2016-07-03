
%{

Debugging notes :

- length drops from ~92k samples for a 9.2s trial to 51k samples for a 7.5
second trial at trial 147.
- no big discontinuities in the signal suggest short stopping.

- in eg. trial 146, there's a ~5ms delay between bpod's recording of punish
  and nidaq's report via Ai1/Ai2 (that said, this time is diff'd so off by
  0.5ms)


Possible (shitty) solutions:
 1 - Downsample to 1khz or 5khz instead of 10khz.  We're getting nothing useful out of 10.
 2 - Try to debug the data ready callback, since that's clearly being weird.
 

NOISE NOTES

- There's significant noise at 34hz, 88hz, and a few others.


%}

%% MISC DEBUGGING CODE




%%


%file_directory = '/Volumes/Macintosh HD/Users/avaughan/Dropbox/KepecsLab (1)/_Alex/Photometry/bpod data/';
file_directory = '~/Documents/MATLAB/_projects_CSHL/photometry_analysis/bpod_data/';

file_names = {};
trials_to_skip = [];
trial_outcomes = {};
notch_frequency = NaN;


%% Early training trials.


%file_names{end+1} = 'ACh_wt1_AudGoNoGo_NIDAQ_Jul14_2015_Session1.mat';

%file_names{end+1} = 'VIPGCaMP1_AudGoNoGo_NIDAQ_Jul14_2015_Session7.mat';


%file_names{end+1} = 'VIPGCaMP1_AudGoNoGo_NIDAQ_Jul15_2015_Session2.mat';
%trials_to_skip = [22];

%file_names{end+1} = 'ACh_wt2_AudGoNoGo_NIDAQ_Jul15_2015_Session3.mat';

% NO LICKING
%file_names{end+1} = 'ACh_wt1_AudGoNoGo_NIDAQ_Jul16_2015_Session1.mat'

% 50 licks %
%file_names{end+1} = 'ACh_wt2_AudGoNoGo_NIDAQ_Jul16_2015_Session1.mat'

% 15 licks %
% file_names{end+1} = 'ACh_wt2_AudGoNoGo_NIDAQ_Jul16_2015_Session2.mat'

% WOOOO
%file_names{end+1} = 'VIPGCaMP1_AudGoNoGo_NIDAQ_Jul16_2015_Session3.mat';
%file_names{end+1} = 'VIPGCaMP2_AudGoNoGo_NIDAQ_Jul16_2015_Session2.mat'


% file_names{end+1} = 'VIPGCaMP2_AudGoNoGo_NIDAQ_Jul19_2015_Session1.mat';
%file_names{end+1} = 'VIPGCaMP1_AudGoNoGo_NIDAQ_Jul19_2015_Session1.mat';
%file_names{end+1} = 'ACh_wt2_AudGoNoGo_NIDAQ_Jul17_2015_Session1.mat';
%file_names{end+1} = 'ACh_wt1_AudGoNoGo_NIDAQ_Jul17_2015_Session3.mat';
% file_names{end+1} = 'VIPGCaMP1_AudGoNoGo_NIDAQ_Jul17_2015_Session3.mat';
% 
% trial_legends   = { 'Punishment','Reward','No Lick'};
% trial_outcomes = { [0],[1],[3] };
% trial_outcome_legends   = { 'Punishment','Reward','Omission'};
% trial_colors    = {'r','g','b','k','k','k'};
% trial_outcome_fields = {'DeliverPunish','DeliverReward','PostTrialRecording'};

%%%%%% %%%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%%
%%%%%% %%%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%%
%%%%%% %%%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%%
%%%%%% %%%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%%
%%%%%% %%%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%%



%% %% %% %% %% %% %% Pavlovian!% %% %% %% %% %% %% %% %% %% %%

%%file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Jul27_2015_Session3.mat';

% July 28
%file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Jul28_2015_Session1.mat';%
%file_names{end+1} = 'VIP_Ach_black_SO_Training_NIDAQ_Jul28_2015_Session2.mat';
%file_names{end+1} = 'VIP_Ach_Pink_SO_Training_NIDAQ_Jul28_2015_Session3.mat';
%file_names{end+1}   = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Jul28_2015_Session1.mat';


% July 29
% This file format is fucked for some reason - try redownloading?
% NO TRIALS % file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Jul29_2015_Session1.mat';

%file_names{end+1} = 'ACh_wt1_black_SO_Training_NIDAQ_Jul29_2015_Session1.mat';

% WEIRD TRIALS % file_names{end+1} = 'VIP_WT_pink_SO_Training_NIDAQ_Jul29_2015_Session1.mat';
%file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Jul29_2015_Session2.mat';

% OK, but punishment not on %
%file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Jul29_2015_Session1.mat';

% No punishment - I don't see a signal either.
%file_names{end+1} = 'VIP_Ach_Pink_SO_Training_NIDAQ_Jul29_2015_Session1.mat';


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% WITH PUNISHMENT %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%

% FEW TRIALS - otherwise ok but no signal.
%file_names{end+1} = 'VIP_Ach_black_SO_Training_NIDAQ_Jul29_2015_Session1.mat';
% VERY FEW TRIALS - otherwise ok but no signal.
%file_names{end+1} = 'VIP_Ach_Pink_SO_Training_NIDAQ_Jul29_2015_Session2.mat';

% Jul 30
% Lots of cutoffs %
% Nice data, though - unexpected reward >> expected reward.
%file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Jul30_2015_Session1.mat'; % NO AUDIO SIGNAL (but probably not hooked up)

% Jul 30  - OK
%file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Jul30_2015_Session2.mat';

% July 30 ACH
% Ok, but no signal.
%file_names{end+1} = 'VIP_Ach_black_SO_Training_NIDAQ_Jul30_2015_Session2.mat';

% FEW TRIALS / No obvious signal
%file_names{end+1} = 'VIP_Ach_Pink_SO_Training_NIDAQ_Jul30_2015_Session2.mat';
% FEW TRIALS / No obvious signal
%file_names{end+1} = 'VIP_Ach_Pink_SO_Training_NIDAQ_Jul30_2015_Session3.mat';

%%%%% Hyun's week of training - not sure any of these are useful??
%%% - many drops
%%% file_names{end+1} = 'VIP_Ach_black_SO_Training_NIDAQ_Aug14_2015_Session1.mat';
%%% - many drops
%%% - seems to crash things on trial 71 or so.
%file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Aug11_2015_Session1.mat';
%%% - Lockin turned OFF
%%% - OK but only 3 trial types?  Also, probably needs MUCH stronger filtering.
% file_names{end+1} = 'VIP_Ach_Pink_SO_Training_NIDAQ_Aug14_2015_Session1.mat';
% file_names{end+1} = 'VIP_Ach_Pink_SO_Training_NIDAQ_Aug12_2015_Session2.mat';

%%% Restarting training Aug 18
%%% - Many drops after trial 151 and missing outcome 1,3,6.
%%% - Lockin not normalized.
%%% - Filtering at 3ms, which makes it much noisier.
%file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Aug18_2015_Session13.mat'; % NO AUDIO SIGNAL

%%% Aug 19 GCaMP 1.
%%% Very noisy - only had 3ms time constant, I think.
%%% Works well enough on lowpass, though with Fpass = 10; Fstop = 20; Ap = 1; Ast = 30;
%file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Aug19_2015_Session1.mat';

%%% Aug 19 GCaMP_2
%%% Set with 10ms time constant for smoothing.
%%% Looks fine with lowpass Fpass = 10; Fstop = 20; Ap = 1; Ast = 30;
%%% But might be overfiltering at that point.
%file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Aug19_2015_Session1.mat';

% Aug 20 -
%file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Aug20_2015_Session2.mat';   % NO AUDIO SIGNAL
%file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Aug20_2015_Session1.mat';       % NO AUDIO SIGNAL

% Aug 21 -
% Most trials today had chopper OFF so are ~useless.  Toss everything but Ach_pink
%file_names{end+1} = 'VIP_Ach_Pink_SO_Training_NIDAQ_Aug21_2015_Session1.mat';      % NO AUDIO SIGNAL


% Aug 26-Sep 01
% Audio may be unplugged in these trials, but they should definitely show 
% audio response in Ai2 (ie., channel 3)
%file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Aug21_2015_Session1.mat';   % NO AUDIO SIGNAL
%file_names{end+1} = 'VIP_Ach_black_SO_Training_NIDAQ_Aug26_2015_Session1.mat';     % NO AUDIO SIGNAL
%file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Aug26_2015_Session1.mat';   % NO AUDIO SIGNAL
%file_names{end+1} = 'VIP_Ach_black_SO_Training_NIDAQ_Aug27_2015_Session1.mat';     % NO AUDIO SIGNAL
%file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Aug27_2015_Session1.mat';   % NO AUDIO SIGNAL
%file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Aug28_2015_Session1.mat';   % NO AUDIO SIGNAL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Sep 01
% At some point before this, audio was unplugged :|
% However, we can see audio in channel 3 for confirmation.

% GCaMP2 (pink) - needs 20hz filter
%  file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Sep01_2015_Session2.mat';
%  file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Sep02_2015_Session1.mat';
 file_names{end+1} = 'Chat_GCaMP_pink_SO_Training_NIDAQ_Sep18_2015_Session1.mat';
% file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Sep03_2015_Session1.mat'; % very bad noise
% file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Sep10_2015_Session2.mat'; % very bad noise
%notch_frequency = 20;

% GCaMP1 (black) - needs
% WHERE IS SEP 01?
% file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Sep02_2015_Session1.mat';
%file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Sep03_2015_Session1.mat';  % Very noisy
%file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Sep10_2015_Session1.mat';  % --> Very weird.
% file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Sep11_2015_Session1.mat';
notch_frequency = 13; % Some 13hz for Sep02/03, but stronger 20hz for Sep10.


%file_names{end+1} = 'VIP_Ach_black_SO_Training_NIDAQ_Sep01_2015_Session1.mat';     % AUDIO OK

%%%%%% %%%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%%
%%%%%% %%%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%%
%%%%%% %%%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%%
%%%%%% %%%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%%
%%%%%% %%%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%%


%% Trial groupings for trials with punishment - 
% *** FOR SO TASKONLY ***

% Trial Groupings are based on actual outcome (reward/punishment) and
% pre-lick or no-pre-lick.

% 0 - pre-lick + PUNISH
% 4 -  no-lick + PUNISH
%
% 1 - pre-lick + REWARD
% 2 -  no-lick + REWARD
%
% 5 - pre-lick + OMISSION
% 3 -  no-lick + OMISSION

% With trial types, these arrive in teh following order
% E(Reward) > Punish
% E(Punish) > Punish
% E(Reward) > Reward
% E(Punish) > Reward
% E(Reward) > Omission
% E(Punish) > Omission
 
trial_outcomes        = {  [0 4],[1 2], [3 5]};
trial_outcome_fields  = { 'Punish',     'Reward', 'PostUS'      };
trial_outcome_legends = { 'Punishment', 'Reward', 'Omission'    };
trial_colors          = { 'm', 'r' ,  'g','c',      'b','k'           };
%%%trial_colors          = { 'r',          'g',      'b'           };

% trial_outcomes = { 0,4,1,2,3,5 };
% trial_feedback_fields = { 'Punish', 'Punish',    'Reward','Reward', 'PostUS', 'PostUS'      };
% trial_legends         = { 'Punish (pre-lick)', 'Punish (no-lick)',    'Reward (pre-lick)','Reward (no-lick)', 'Omission (pre-lick)', 'Omission (no-lick)'      };
% trial_colors          = { 'm', 'r' ,  'c','g',      'b','k'           };

% Trial types are 1 (predominantly reward) and 0 (predominantly punishment)
trial_types        = [1 2]; % <-- drawn from data below, but this is usually what
trial_type_legends = {'E(Reward)','E(Punish)'};
%it is.


%% Analysis parameters

% DO ANALYSIS OR JUST PRETENDE
actually_do_analysis            = 1;

% COMBINE ALL FILES ON A TRIAL-BY-TRIAL BASIS?
combine_trials_across_sessions  = 0;

% Do all plots, or only sanity-check plots?
do_all_plots = 1;

% Do all pairwise plots?  This can get annoying.
do_all_pairwise_plots           = 0;


% SKIP SOME TRIALS FOR ALL FILES?
% Only Trials 5:55
%trials_to_skip_for_all_fn = @(n_trials) [1:5 56:n_trials n_trials-[0:4]];
% Skip 1st 5 and last 5 trials, and anything over trial 250
trials_to_skip_for_all_fn = @(n_trials) [1:5, n_trials-[0:4], 251:1000];

% Filtering?
do_filtering = 0;
filter_type = {'lowpass','notch'}; % 'notch' or 'lowpass', 'bandpass'

% Time (s) on the either side of the alignment cue
% _range is the full range for filtering, etc.
% _display is the actual xlims on display
psth_flanking_range   = [-2.0 2.1]; % max pre is ~2s because of stimulus presentation.  Can be changed though
psth_flanking_xlim   = [-1.5 2];


%%
assert(~isempty(file_names),'file_names is empty :: You should probably select some files to analyze')


sample_rate = 10000;
% Size of window for zeroing response.
zeroing_samples = sample_rate;

% FILTERING
% Low pass
Fpass = 45; Fstop = 50; Ap = 2; Ast = 20;
d_lowpass = designfilt('lowpassiir','PassbandFrequency',Fpass,...
    'StopbandFrequency',Fstop,'PassbandRipple',Ap,...
    'StopbandAttenuation',Ast,'SampleRate',sample_rate);

% Notch (bandstop)
if isnan(notch_frequency)
    notch_frequency = 20; % or 34...
end
d_notch = designfilt('bandstopiir',...
    'FilterOrder',20, ...
    'HalfPowerFrequency1',notch_frequency * 0.9,...
    'HalfPowerFrequency2',notch_frequency * 1.1, ...
    'SampleRate',sample_rate);

% d_notch = designfilt('bandstopiir','FilterOrder',20, ...
%     'HalfPowerFrequency1',35,'HalfPowerFrequency2',45, ...
%     'SampleRate',sample_rate);

% Bandpass.
d_bandpass = designfilt('bandpassiir',...
    'FilterOrder',20, ...
    'HalfPowerFrequency1',0.5,...
    'HalfPowerFrequency2',30, ...
    'SampleRate',sample_rate);

%%
switch do_filtering
    case 1
        filter_string = ' :: FILTERING ON';
    case 0
        filter_string = ' :: FILTERING OFF';
end


% Y range for all figures
ylims_all = [-2 2];


for file_number = 1:length(file_names);
    
    % Load file
    file_name = file_names{file_number};
    load([file_directory file_name]);
    
    fprintf('\nLoaded file %s\n',file_name)
    
    % Trial params
    n_trials = length(SessionData.TrialTypes);
    n_channels = size(SessionData.NidaqData{1},2);
    fluorescent_channel = 0 +1; % For Ai0
    audio_channel = [2 3];
    
    % PSTH samples.
    psth_flanking_times   = linspace( psth_flanking_range(1),psth_flanking_range(2), sample_rate*diff(psth_flanking_range)+1  );
    psth_flanking_samples = psth_flanking_times * sample_rate;
    
    % Default if not specified above.
    if isempty(trial_outcomes)
        trial_outcomes = num2cell(unique(SessionData.TrialOutcome));
    end
    
    %% Generate nd_trial_typees and nd_trial_outcomes to iterate over.
    trial_types    = unique(SessionData.TrialTypes);
    
    % Index
    [nd_trial_types,nd_outcome_types] = ndgrid(1:length(trial_types),1:length(trial_outcomes));
    trial_type_outcome_inds = [nd_trial_types(:) nd_outcome_types(:)];
    n_trial_type_vs_outcomes = size(trial_type_outcome_inds,1);
    
    %%
    fprintf('Available trial types: %s\n', mat2str(cell2mat(trial_outcomes)))
    fprintf('Available trial outcomes: %s\n',mat2str(cell2mat(trial_outcomes)))
    fprintf(' %.0f pairwise combinations\n',n_trial_type_vs_outcomes)
    
    %% Trial Size Figure and cutoff for short trials
    namedFigure(sprintf('Trial Sizes :: %s',file_name)); clf; hold on
    [trial_sizes, postUS_times] = deal([]);
    for i = 1:SessionData.nTrials, trial_sizes(i) = size(SessionData.NidaqData{i},1)/sample_rate;,end
    for i = 1:SessionData.nTrials,
        try,
            % SO (stimulus-outcome_ training
            postUS_times(i) = SessionData.RawEvents.Trial{i}.States.PostUS(1);
        catch
            % AudGoNoGo task.
            postUS_times(i) = SessionData.RawEvents.Trial{i}.States.WaitForLick(2);
        end
    end
    plot(trial_sizes,'b')
    plot(postUS_times,'g')
    plot(postUS_times > trial_sizes,'r.')
    legend({'Recording Length','PostUS Time','MISMATCH'})
    
    trial_too_short_cutoff = find(postUS_times > trial_sizes,1);
    if isempty(trial_too_short_cutoff), trial_too_short_cutoff = SessionData.nTrials+1; end
    fprintf('trial_too_short_cutoff is set to trial #%.0f\n',trial_too_short_cutoff)
    
    %% Bleaching curve.
    namedFigure(sprintf('Bleaching Curve :: %s',file_name)); clf
    subplot(2,2,1); cla; hold on
    for i = 1:SessionData.nTrials
        start_times(i) = (SessionData.TrialStartTimestamp(i) - SessionData.TrialStartTimestamp(1))/60;
        start_fluor(i) = mean(SessionData.NidaqData{i}(1:sample_rate));
    end
    h1 = plot( start_times,start_fluor,'o');
    xlims = xlim;
    x_linspace = linspace(xlims(1),xlims(2),100);
    h2 = plot( x_linspace, start_fluor(1) * 0.99.^(x_linspace),'g:');
    h3 = plot( x_linspace, start_fluor(1) * 0.98.^(x_linspace),'r:');
    legend([h1,h2,h3],{'Observed','1%/min bleaching rate','2%/min bleaching rate'})
    xlabel('Time (minutes?)')
    ylabel('Raw signal')
    title('Raw fluorescence in first 1s of trial')
    subplot(2,2,2); cla; hold on
    for i = 1:SessionData.nTrials
        start_times(i) = (SessionData.TrialStartTimestamp(i) - SessionData.TrialStartTimestamp(1))/60;
        start_fluor(i) = mean(SessionData.NidaqData{i}(1:sample_rate)) / mean(SessionData.NidaqData{1}(1:sample_rate));
    end
    h1 = plot( start_times,start_fluor,'o');
    xlims = xlim;
    x_linspace = linspace(xlims(1),xlims(2),100);
    h2 = plot( x_linspace, 0.99.^(x_linspace),'g:');
    h3 = plot( x_linspace, 0.98.^(x_linspace),'r:');
    legend([h1,h2,h3],{'Observed','1%/min bleaching rate','2%/min bleaching rate'})
    xlabel('Time (minutes?)')
    ylabel('Normalized signal')
    title({'Note :: bleaching rate not valid','if lockin amplifier was pre-zeroed.'})
    subplot(2,2,3); cla; hold on
    for i = 1:SessionData.nTrials
        start_fluor(i) = mean(SessionData.NidaqData{i}(1:sample_rate));
    end
    h1 = plot(start_fluor,'o');
    xlabel('Trials')
    ylabel('Raw signal')
    title('Raw fluorescence in first 1s of trial')
    subplot(2,2,4); cla; hold on
    for i = 1:SessionData.nTrials
        start_fluor(i) = mean(SessionData.NidaqData{i}(1:sample_rate)) / mean(SessionData.NidaqData{1}(1:sample_rate));
    end
    h1 = plot( start_fluor,'o');
    xlabel('Trials')
    ylabel('Normalized signal c')
    title({'Note :: bleaching rate not valid','if lockin amplifier was pre-zeroed.'})
    
    %%
    
    if actually_do_analysis
        
        if combine_trials_across_sessions
            audio_fig                   = namedFigure(['Combined__audio' filter_string]); clf
            if do_all_plots
                average_fig                 = namedFigure(['Combined__average' filter_string]); clf
                all_trials_fig              = namedFigure(['Combined__all_trials' filter_string]); clf
                spectrogram_fig             = namedFigure(['Combined__spectrogram' filter_string]); clf
                combined_fig                = namedFigure(['Combined__combined' filter_string]); clf
                sample_fig                  = namedFigure(['Combined__sample' filter_string]); clf
                comparison_full_fig         = namedFigure(['Combined__comparison_figure' filter_string]); clf
                lick_raster_fig             = namedFigure(['Combined__lick_raster' filter_string]); clf
                if do_all_pairwise_plots
                    pairwise_stim_fig       = namedFigure(['Combined__pairwise_stim' filter_string]); clf
                    pairwise_feedback_fig   = namedFigure(['Combined__pairwise_feedback' filter_string]); clf
                end
            end
        else
            audio_fig                   = namedFigure([file_name '__audio' filter_string]); clf
            if do_all_plots
                average_fig                 = namedFigure([file_name '__average' filter_string]); clf
                all_trials_fig              = namedFigure([file_name '__all_trials' filter_string]); clf
                spectrogram_fig             = namedFigure([file_name '__spectrogram' filter_string]); clf
                combined_fig                = namedFigure([file_name '__combined' filter_string]); clf
                sample_fig                  = namedFigure([file_name '__sample' filter_string]); clf
                pairwise_stim_fig           = namedFigure([file_name '__pairwise_stim' filter_string]); clf
                pairwise_feedback_fig       = namedFigure([file_name '__pairwise_feedback' filter_string]); clf
                comparison_full_fig         = namedFigure([file_name '__comparison_figure' filter_string]); clf
                lick_raster_fig             = namedFigure([file_name '__lick_raster' filter_string]); clf
                if do_all_pairwise_plots
                    pairwise_stim_fig       = namedFigure([file_name '__pairwise_stim' filter_string]); clf
                    pairwise_feedback_fig   = namedFigure([file_name '__pairwise_feedback' filter_string]); clf
                end
            end
            
        end
        
    end
    
    
    % Loop across trial types
    sample_trials = cell(size(n_trial_type_vs_outcomes));
    for trial_type_outcome_ind  = 1:length(trial_type_outcome_inds ),
        
        %%
        trial_type_ind = trial_type_outcome_inds(trial_type_outcome_ind,1);
        trial_outcome_ind = trial_type_outcome_inds(trial_type_outcome_ind,2);
        trial_type    = trial_types(trial_type_ind);
        trial_outcome = trial_outcomes{trial_outcome_ind};
        
        trial_type_outcome_legend = sprintf('%s > %s',trial_type_legends{trial_type_ind},trial_outcome_legends{trial_outcome_ind});
        
        fprintf('trial_type %.0f && trial_outcome %s\n',trial_type,mat2str(trial_outcome))
        
        if length(trial_outcome) > 1
            % Group more than one condition into the same bin
            trial_group_indices = find(sum( bsxfun(@eq,SessionData.TrialOutcome,trial_outcome')));
        else
            % Find onyl a single condition
            trial_group_indices = find(bsxfun(@eq,SessionData.TrialOutcome,trial_outcome'));
        end
        fprintf('  trial_outcome = %s \n',mat2str(trial_outcome))
        fprintf('    %03.0f trials\n',length(trial_group_indices))
        
        fprintf('  culling by trial_type = %.0f\n',trial_type)
        trial_group_indices = intersect(trial_group_indices,find(SessionData.TrialTypes == trial_type));
        fprintf('    %03.0f trials remaining\n',length(trial_group_indices))
        
        fprintf('  culling by trials dropped that exceed trial_too_short_cutoff\n')
        trial_group_indices  = trial_group_indices(trial_group_indices < trial_too_short_cutoff);
        fprintf('    %03.0f trials remaining .\n',length(trial_group_indices))
        
        fprintf('  trimming with trials_to_skip and trials_to_skip_for_all_fn\n')
        trial_group_indices = setdiff(trial_group_indices,trials_to_skip);
        trial_group_indices = setdiff(trial_group_indices,trials_to_skip_for_all_fn(n_trials));
        trial_group_indices = round(trial_group_indices);
        fprintf('    %03.0f trials remaining.\n',length(trial_group_indices))
        
        if length( trial_group_indices ) < 2
            continue
        end
        
        if ~actually_do_analysis
            continue
        end
        
        %%
        % This determines which ROW of the assembeld response matrix we will draw from
        % -- NOT which trial number it is.
        %sample_trials{trial_group_ind} = ceil(rand * numel(trial_group_indices));
        
        % Loop through trials for this trial type - we do it this way to
        % allow expansion when combining across different sessions.
        trial_num = 1;
        
        % TODO :: FIX THIS - 
        if not(combine_trials_across_sessions) || file_number == 1
            [peri_stim_GCaMP{trial_type_outcome_ind}, peri_feedback_GCaMP{trial_type_outcome_ind}] = deal(     nans(length(trial_group_indices),length(psth_flanking_samples))   );
            [peri_stim_audio{trial_type_outcome_ind}, peri_feedback_audio{trial_type_outcome_ind}] = deal(     nans(length(trial_group_indices),length(psth_flanking_samples))   );
            [peri_stim_licks{trial_type_outcome_ind}, peri_feedback_licks{trial_type_outcome_ind}] = deal(   sparse(length(trial_group_indices),length(psth_flanking_samples))  );
            
        end
        
        
        for trial = trial_group_indices,
            
            %% Stimulus
            
            % GCaMP response
            stimulus_start_time = SessionData.RawEvents.Trial{trial}.States.DeliverStimulus(1);
            sample_indices = round( stimulus_start_time * sample_rate + psth_flanking_samples(:) );
            available_samples = size( SessionData.NidaqData{trial} , 1) - stimulus_start_time * sample_rate;
            if available_samples < 0,
                fprintf('available_samples < 0 for some reason (%.0f) :: trial %.0f :: stimulus\n',available_samples,trial)
                continue
            end
            end_sample = round(min(length(psth_flanking_samples),available_samples));
            % Symmetric window around start
            %stimulus_zeroing_window = (stimulus_start_time * sample_rate) + round(linspace(-floor(zeroing_samples/2),floor(zeroing_samples/2)-1,zeroing_samples));
            % Causal window around start
            stimulus_zeroing_window = (stimulus_start_time * sample_rate) + round(linspace(-floor(zeroing_samples),-1,zeroing_samples));
            % Populate peri_stim_GCaMP matrix with response normalized by mean around t=0;
            peri_stim_GCaMP{trial_type_outcome_ind}( trial_num, 1:end_sample ) = ...
                SessionData.NidaqData{trial}( sample_indices(1:end_sample) , fluorescent_channel )' ...
                - mean(SessionData.NidaqData{trial}( round(stimulus_zeroing_window) , fluorescent_channel ));
            
            peri_stim_audio{trial_type_outcome_ind}( trial_num, 1:end_sample ) = ...
                mean(SessionData.NidaqData{trial}( sample_indices(1:end_sample) , audio_channel ),2)' ...
                - mean(mean(SessionData.NidaqData{trial}( round(stimulus_zeroing_window) , audio_channel),2));
            
            
            
            %% Feedback GCaMP Response
            if trial_outcome_ind> length(trial_outcome_fields)
                fprintf('ERROR :: trial_outcome_ind > length(trial_outcome_fields)\n')
                keyboard
                continue
                
            else
                
                feedback_start_time = SessionData.RawEvents.Trial{trial}.States.(trial_outcome_fields{trial_outcome_ind})(1);
                
                % If ~punish --> no punish
                % If we don't think trial is a punish trial, then there should be no punishment!
                % State 'Punish' is from the SO_ task, 'DeliverPunish' is from AudGoNoGo
                if isfield(SessionData.RawEvents.Trial{trial}.States,'Punish')
                    punish_field = SessionData.RawEvents.Trial{trial}.States.Punish;
                elseif isfield(SessionData.RawEvents.Trial{trial}.States,'DeliverPunish')
                    punish_field = SessionData.RawEvents.Trial{trial}.States.DeliverPunish;
                else
                    error('Could not find punishment field.')
                end
                assert( xor( sum(strcmp(trial_outcome_fields{trial_outcome_ind}, {'Punish','DeliverPunish'})) , ...
                    isnan( punish_field(1))   ));
                
                % If ~reward trial --> no reward
                % If we don't think trial is a reward trial, then there should be no reward!
                % State 'Reward' is from the SO_ task, 'DeliverReward' is from AudGoNoGo
                if isfield(SessionData.RawEvents.Trial{trial}.States,'Reward')
                    reward_field = SessionData.RawEvents.Trial{trial}.States.Reward;
                elseif isfield(SessionData.RawEvents.Trial{trial}.States,'DeliverReward')
                    reward_field = SessionData.RawEvents.Trial{trial}.States.DeliverReward;
                else
                    error('Could not find Reward field.')
                end
                assert( xor( sum(strcmp(trial_outcome_fields{trial_outcome_ind}, {'Reward','DeliverReward'})) , ...
                    isnan( reward_field(1))   ));
                
                
            end
            assert( ~isnan(feedback_start_time) , 'Feedback start time is a NaN - usually a sign that you''re looking at the wrong behavioral timepoint.')
            sample_indices = round( feedback_start_time * sample_rate + psth_flanking_samples(:) );
            available_samples = size( SessionData.NidaqData{trial} , 1) - feedback_start_time * sample_rate;
            if available_samples < 0,
                fprintf('available_samples < 0 for some reason :: trial %.0f :: reward\n',trial)
                continue
            end
            
            % This will catch trials that ended too early.
            %             if (available_samples / sample_rate) < 2,
            %                 keyboard
            %             end
            
            end_sample = round(min(length(psth_flanking_samples),available_samples));
            
            % Normalizing by feedback start time to be dF/F is not very
            % helpful the baseline value is below zero.  So we just do dF.
            feedback_zeroing_window = feedback_start_time * sample_rate + linspace(-floor(zeroing_samples/2),floor(zeroing_samples/2));
            peri_feedback_GCaMP{trial_type_outcome_ind}( trial_num, 1:end_sample ) = ...
                ( SessionData.NidaqData{trial}( sample_indices(1:end_sample) , fluorescent_channel )' ...
                - mean(SessionData.NidaqData{trial}( round(feedback_zeroing_window) , fluorescent_channel)));
            
            peri_feedback_audio{trial_type_outcome_ind}( trial_num, 1:end_sample ) = ...
                mean( SessionData.NidaqData{trial}( sample_indices(1:end_sample) , audio_channel ),2)' ...
                - mean(mean(SessionData.NidaqData{trial}( round(feedback_zeroing_window) , audio_channel),2));
            %/ SessionData.NidaqData{trial}( round(feedback_start_time * sample_rate) , fluorescent_channel ) ;
            
            
            %% Filtering
            
            % Notch or smoothing filtering.
            % Make sure we only filter the parts of this matrix that are
            % sensible - outside of 1:end_sample may be NaNs at this point.
            peri_stim_GCaMP_filt{trial_type_outcome_ind} = peri_stim_GCaMP{trial_type_outcome_ind};
            if do_filtering
                
                if any(strcmp(filter_type,'bandpass'))
                    peri_stim_GCaMP_filt{trial_type_outcome_ind}( trial_num, 1:end_sample )  = filtfilt(d_bandpass, peri_stim_GCaMP_filt{trial_type_outcome_ind}(trial_num, 1:end_sample)');
                end
                if any(strcmp(filter_type,'lowpass'))
                    peri_stim_GCaMP_filt{trial_type_outcome_ind}( trial_num, 1:end_sample )  = filtfilt(d_lowpass, peri_stim_GCaMP_filt{trial_type_outcome_ind}(trial_num, 1:end_sample)');
                end
                if any(strcmp(filter_type,'notch'))
                    peri_stim_GCaMP_filt{trial_type_outcome_ind}( trial_num, 1:end_sample )  = filtfilt(d_notch, peri_stim_GCaMP_filt{trial_type_outcome_ind}(trial_num, 1:end_sample)');
                end
            end
            
            % Notch or smoothing filtering.
            % Make sure we only filter the parts of this matrix that are
            % sensible - outside of 1:end_sample may be NaNs at this point.
            peri_feedback_GCaMP_filt{trial_type_outcome_ind} = peri_feedback_GCaMP{trial_type_outcome_ind};
            if do_filtering
                if any(strcmp(filter_type,'bandpass'))
                    peri_feedback_GCaMP_filt{trial_type_outcome_ind}( trial_num, 1:end_sample )  = filtfilt(d_bandpass, peri_feedback_GCaMP_filt{trial_type_outcome_ind}(trial_num, 1:end_sample)');
                end
                if any(strcmp(filter_type,'lowpass'))
                    peri_feedback_GCaMP_filt{trial_type_outcome_ind}( trial_num, 1:end_sample )  = filtfilt(d_lowpass, peri_feedback_GCaMP_filt{trial_type_outcome_ind}(trial_num, 1:end_sample)');
                end
                if any(strcmp(filter_type,'notch'))
                    peri_feedback_GCaMP_filt{trial_type_outcome_ind}( trial_num, 1:end_sample )  = filtfilt(d_notch, peri_feedback_GCaMP_filt{trial_type_outcome_ind}(trial_num, 1:end_sample)');
                end
            end
            
            
            %% Calculate lick rate.
            
            % Stimulus
            stimulus_start_time = SessionData.RawEvents.Trial{trial}.States.DeliverStimulus(1);
            sample_indices = round( stimulus_start_time * sample_rate + psth_flanking_samples(:) );
            available_samples = size( SessionData.NidaqData{trial} , 1) - stimulus_start_time * sample_rate;
            
            % Gather licks.
            if isfield(SessionData.RawEvents.Trial{trial}.Events,'Port2In')
                all_licks = SessionData.RawEvents.Trial{trial}.Events.Port2In;
                fprintf('    Trial %.0f :: %.0f Port2In events\n',trial,length(all_licks ))
            else
                fprintf('    Trial %.0f :: No Port2In events!\n',trial)
                all_licks = [];
            end
            
            % Align to PSTH
            stimulus_window = stimulus_start_time + psth_flanking_range;
            feedback_window = feedback_start_time + psth_flanking_range;
            
            this_peri_stim_licks     = all_licks( all_licks>stimulus_window(1) & all_licks < stimulus_window(2));
            this_peri_feedback_licks = all_licks( all_licks>feedback_window (1) & all_licks < feedback_window (2));
            
            % Convert to samples after beginning of PSTH
            this_peri_stim_licks_psth_samples = ceil(sample_rate*(this_peri_stim_licks-(stimulus_start_time+psth_flanking_range(1))));
            this_peri_feedback_licks_psth_samples = ceil(sample_rate*(this_peri_feedback_licks-(feedback_start_time+psth_flanking_range(1))));
            
            % Make a big (sparse!) matrix for this.
            peri_stim_licks{trial_type_outcome_ind}(trial_num, this_peri_stim_licks_psth_samples ) = 1;
            peri_feedback_licks{trial_type_outcome_ind}(trial_num,this_peri_feedback_licks_psth_samples) = 1;
            
            
            %% Advance trial_num
            trial_num = trial_num + 1;
            
        end % Iterating over trial_group index
        
        
        %% Error handling on the length of the peri_feedback_GCaMP etc.
        if size(peri_feedback_GCaMP_filt{trial_type_outcome_ind},2)  < length(psth_flanking_samples)
            fprintf('  *** WARNING *** :: trial_type_outcome_ind %.0f :: expected %.0f samples but only found %.0f\np',...
                trial_type_outcome_ind,...
                length(psth_flanking_samples),...
                size(peri_feedback_GCaMP_filt{trial_type_outcome_ind},2))
        end
        
        
        
        
        %% If skipping plotting, continue looping here.
        
        if combine_trials_across_sessions && ~(file_number == length(file_names))
            fprintf('Skipping plots for now - combining later.\n')
            continue
        end
        
        
        
        %% Prep for plotting
        
        % Find representative sample (ie., most like the average);
        %peri_feedback_GCaMP_filt_nansafe = peri_feedback_GCaMP_filt{trial_group_ind};
        %peri_feedback_GCaMP_filt_nansafe(isnan(peri_feedback_GCaMP_filt_nansafe )) = 0;
        template_errors = sum(bsxfun(@minus,peri_feedback_GCaMP_filt{trial_type_outcome_ind},nanmean(peri_feedback_GCaMP_filt{trial_type_outcome_ind})).^2,2);
        sample_trials{trial_type_outcome_ind} = find(template_errors == min(template_errors),1);
        
        
        % Unfiltered - includes all trials - nansafe
        peri_stim_GCaMP_to_plot     =     peri_stim_GCaMP{trial_type_outcome_ind};
        peri_feedback_GCaMP_to_plot = peri_feedback_GCaMP{trial_type_outcome_ind};
        
        peri_stim_GCaMP_to_plot(isnan(peri_stim_GCaMP_to_plot))         = 0;
        peri_feedback_GCaMP_to_plot(isnan(peri_feedback_GCaMP_to_plot)) = 0;
        
        %Filtered - includes all trials - nansafe
        peri_stim_GCaMP_filt_to_plot     =     peri_stim_GCaMP_filt{trial_type_outcome_ind};
        peri_feedback_GCaMP_filt_to_plot = peri_feedback_GCaMP_filt{trial_type_outcome_ind};
        
        peri_stim_GCaMP_filt_to_plot(isnan(peri_stim_GCaMP_filt_to_plot)) = 0;
        peri_feedback_GCaMP_filt_to_plot(isnan(peri_feedback_GCaMP_filt_to_plot)) = 0;
        
        % Filtered - AVERAGE AND STD - nansafe
        peri_stim_mean_filt     = nanmean( peri_stim_GCaMP_filt{trial_type_outcome_ind}, 1     );
        peri_stim_err_filt      = nanstd(  peri_stim_GCaMP_filt{trial_type_outcome_ind}, [], 1 ) / sqrt(length(trial_group_indices));
        peri_stim_n             = size(    peri_stim_GCaMP_filt{trial_type_outcome_ind}, 1     );
        
        peri_feedback_mean_filt = nanmean( peri_feedback_GCaMP_filt{trial_type_outcome_ind}, 1     );
        peri_feedback_err_filt  = nanstd(  peri_feedback_GCaMP_filt{trial_type_outcome_ind}, [], 1 ) / sqrt(length(trial_group_indices));
        peri_feedback_n         = size(    peri_feedback_GCaMP_filt{trial_type_outcome_ind}, 1     );
        
        
        
        %%
        
        
        %% Plot :: Audio
        
        figure(audio_fig)
        subplot( n_trial_type_vs_outcomes,2, 2 * (trial_type_outcome_ind-1) + 1)
        imagesc(...
            (1:size(peri_stim_audio{trial_type_outcome_ind},2))/sample_rate + psth_flanking_range(1),...
            1:size(peri_stim_audio{trial_type_outcome_ind},1),...
            peri_stim_audio{trial_type_outcome_ind})
        caxis(ylims_all); colorbar
        if mod(trial_type_outcome_ind,n_trial_type_vs_outcomes) == 1,
            ylabel(trial_type_outcome_legend)
        end
        xlim(psth_flanking_xlim)
        
        subplot( n_trial_type_vs_outcomes, 2, 2 * (trial_type_outcome_ind-1) + 2)
        imagesc(...
            (1:size(peri_feedback_audio{trial_type_outcome_ind},2))/sample_rate + psth_flanking_range(1),...
            1:size(peri_feedback_audio{trial_type_outcome_ind},1),...
            peri_feedback_audio{trial_type_outcome_ind})
        caxis(ylims_all); colorbar
        colormap red_white_blue
        xlim(psth_flanking_xlim)
        
        
        if ~do_all_plots
            continue
        end
       
        
        %% PLOT :: LICKS :: Raster of each trial

        
        %%
        %if ~isempty(lick_times)
        figure(lick_raster_fig)
        
        % Smoothing kernel
        lick_rate_ylim = [0 10];
        sigma= 1000; X = -3*sigma:3*sigma;
        kernel = 1/sqrt(2*pi)*sigma * exp(-0.5*X.^2/(sigma^2));
        %kernel = ones(1,sample_rate/5);
        
        % Stimulus
        
        timepoint_lick_data = {peri_stim_licks{trial_type_outcome_ind},peri_feedback_licks{trial_type_outcome_ind}};
        timepoint_legend = {'STIMULUS','FEEDBACK'}
        for timepoint_ind = 1:2
            
            this_timpoint_data = timepoint_lick_data{timepoint_ind};
            
            subplot( n_trial_type_vs_outcomes,2, 2 * (trial_type_outcome_ind-1) + timepoint_ind); cla; hold on
            [lick_trials,lick_times] = find(this_timpoint_data );
            if ~isempty(lick_trials)
                convolved_lick_rate = zeros(size(this_timpoint_data ));
                for i = 1:size(this_timpoint_data ,1)
                    convolved_lick_rate(i,:) = conv(full(this_timpoint_data (i,:)),kernel,'same');
                end
                convolved_lick_rate = mean(convolved_lick_rate);
                [hAx,hLine1,hLine2] = plotyy( ...
                    lick_times/sample_rate  + psth_flanking_range(1), lick_trials,...
                    psth_flanking_times, convolved_lick_rate*(sample_rate/sum(kernel))  );
                set(hAx,{'ycolor'},{'k';'r'})
                % First axis
                set(hLine1,'LineStyle','.','Color','k');
                ylim([0 trial_num]+[0.5 -0.5]);
                set(hAx(1),'YTick',unique(round(linspace(1,trial_num-1,3))))
                % Second axis
                hAx2 = get(hLine2,'Parent');
                set(hLine2,'LineStyle','-','Color','r','LineWidth',2);
                ylim(hAx(2),lick_rate_ylim )
                set(hAx(2),'YTick',lick_rate_ylim)
                
                switch timepoint_ind
                    case 1
                        ylabel('Trial #');
                    case 2
                        ylabel(hAx2,{'Click','rate (hz)'},'Rotation',-90,'VerticalAlignment','bottom','Color','r')
                end
                
                xlim(hAx(1),psth_flanking_xlim)
                xlim(hAx(2),psth_flanking_xlim)
                
            end
            axis ij
            vline(0)
            title([ timepoint_legend{timepoint_ind} ' :: ' trial_type_outcome_legend])
        end
    
    %end
    
    %% PLOT :: Sample trial
    if ~isempty(sample_trials{trial_type_outcome_ind})
        figure(sample_fig)
        subplot(1,2,1); hold on
        plot( ...
            psth_flanking_times(1:length(peri_stim_GCaMP_filt{trial_type_outcome_ind})), ...
            peri_stim_GCaMP_filt{trial_type_outcome_ind}( sample_trials{trial_type_outcome_ind}, :),...
            trial_colors{trial_type_outcome_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_type_outcome_ind}(1) );
        title({'Aligned to Stimulus',sprintf('Trial %.0f', sample_trials{trial_type_outcome_ind})})
        Outcome_Label = sprintf('%s (%.0f)',trial_type_outcome_legend ,trial_type_outcome_ind);
        ylabel(Outcome_Label)
        ylim(ylims_all*1.5); vline(0,'k')
        xlim(psth_flanking_xlim)
        subplot(1,2,2); hold on
        plot( ...
            psth_flanking_times(1:length(peri_feedback_GCaMP_filt{trial_type_outcome_ind})),...
            peri_feedback_GCaMP_filt{trial_type_outcome_ind}( sample_trials{trial_type_outcome_ind}, :),...
            trial_colors{trial_type_outcome_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_type_outcome_ind}(1) );
        title({'Aligned to Stimulus',sprintf('Trial %.0f', sample_trials{trial_type_outcome_ind})})
        %ylabel(sprintf('Outcome Type %.0f',trial_group_ind))
        ylim(ylims_all*1.5); vline(0,'k')
        xlim(psth_flanking_xlim)
        
    end
    
    %% PLOT Spectrogram of each trial type
    
    periodogram_fn = @(x) periodogram(nanmean(x,1),...
        hamming(size(x,2)),...
        2.^nextpow2(size(x,2)),sample_rate);
    
    figure(spectrogram_fig); 
    subplot( n_trial_type_vs_outcomes,2, 2 * (trial_type_outcome_ind-1) + 1); cla; hold on
    % Plot raw in red.
    [pxx,f] = periodogram_fn(peri_stim_GCaMP_to_plot);
    plot(f*1E-3,10*log10(pxx),'r')
    % Plot filtered in blue.
    periodogram_fn(peri_stim_GCaMP_filt_to_plot);
    set(gca,'XScale','log');
    if mod(trial_type_outcome_ind,n_trial_type_vs_outcomes) == 1,
        h = ylabel(trial_type_outcome_legend);
        set(h,'Rotation',0,'HorizontalAlignment','right','FontSize',14)
    end
    set(gca,'XScale','log')
    set(gca,'XTickLabel',arrayfun(@(x) sprintf('%.0f',x*1000),get(gca,'XTick'),'UniformOutput',false))
    legend({'Raw','Filtered'})
    
    subplot( n_trial_type_vs_outcomes, 2, 2 * (trial_type_outcome_ind-1) + 2); cla; hold on
    hold on
    [pxx,f] = periodogram_fn(peri_feedback_GCaMP_to_plot); % Unfiltered
    plot(f*1E-3,10*log10(pxx),'r')
    periodogram_fn(peri_feedback_GCaMP_filt_to_plot); % Filtered
    set(gca,'XScale','log')
    set(gca,'XTickLabel',arrayfun(@(x) sprintf('%.0f',x*1000),get(gca,'XTick'),'UniformOutput',false))
    ylabel('')
    legend({'Raw','Filtered'})
    
    %% PLOT :: Raster of each trial
    figure(all_trials_fig)
    subplot( n_trial_type_vs_outcomes,2, 2 * (trial_type_outcome_ind-1) + 1)
    imagesc(...
        psth_flanking_range(1) + (1:size(peri_stim_GCaMP_filt_to_plot,2))/sample_rate,...
        1:size(peri_stim_GCaMP_filt_to_plot,1),...
        peri_stim_GCaMP_filt_to_plot); % Filtered!
    caxis(ylims_all); colorbar
    if mod(trial_type_outcome_ind,n_trial_type_vs_outcomes) == 1,
        ylabel(trial_type_outcome_legend)
    end
    xlim(psth_flanking_xlim)
    
    subplot( n_trial_type_vs_outcomes, 2, 2 * (trial_type_outcome_ind-1) + 2)
    imagesc( ...
        psth_flanking_range(1) + (1:size(peri_feedback_GCaMP_filt_to_plot,2))/sample_rate,...
        1:size(peri_feedback_GCaMP_filt_to_plot,1),...
        peri_feedback_GCaMP_filt_to_plot)
    caxis(ylims_all); colorbar
    colormap red_white_blue
    xlim(psth_flanking_xlim)
    
    % PLOT :: Histogram
    %         namedFigure([file_name '__hist' filter_string]);
    %         subplot( n_trial_type_vs_outcomes,2, 2 * (trial_outcomes-1) + 1)
    %         hist(peri_stim_GCaMP_filt(:),100);
    %         subplot( n_trial_type_vs_outcomes, 2, 2 * (trial_outcomes-1) + 2)
    %         hist(peri_feedback_GCaMP_filt(:),100);
    
    %% PLOT :: peri-stimulus / peri-feedback AVERAGE
    figure(average_fig)
    subplot( n_trial_type_vs_outcomes,2, 2 * (trial_type_outcome_ind-1) + 1); hold on
    shadedErrorBar( ...
        psth_flanking_times(1:length(peri_stim_mean_filt)), ...
        peri_stim_mean_filt,      ...
        peri_stim_err_filt,       ...
        {trial_colors{trial_type_outcome_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_type_outcome_ind}(1)} );
    title({'Aligned to Stimulus',sprintf('%.0f trials', length(trial_group_indices) )})
    
    Outcome_Label = sprintf('%s (%.0f)',trial_type_outcome_legend ,trial_type_outcome_ind);
    ylabel(Outcome_Label)
    ylim(ylims_all); vline(0,'k')
    xlim(psth_flanking_xlim)
    
    
    subplot( n_trial_type_vs_outcomes, 2, 2 * (trial_type_outcome_ind-1) + 2); hold on
    shadedErrorBar( ...
        psth_flanking_times(1:length(peri_feedback_mean_filt)),...
        peri_feedback_mean_filt,   ...
        peri_feedback_err_filt,    ...
        {trial_colors{trial_type_outcome_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_type_outcome_ind}(1)} );
    title({'Aligned to Feedback',sprintf('%.0f trials', length(trial_group_indices) )})
    %ylabel(sprintf('Outcome Type %.0f',trial_group_ind))
    ylim(ylims_all); vline(0,'k')
    xlim(psth_flanking_xlim)
    
    
    %% PLOT :: peri-stimulus average combined figure
    
    figure(combined_fig)
    subplot( 1, 2, 1); hold on;
    shadedErrorBar( ...
        psth_flanking_times(1:length(peri_stim_mean_filt)), ...
        peri_stim_mean_filt,      ...
        peri_stim_err_filt,       ...
        {trial_colors{trial_type_outcome_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_type_outcome_ind}(1)} );
    title({'Aligned to Stimulus'})
    Outcome_Label = sprintf('%s (%.0f)',trial_type_outcome_legend ,trial_type_outcome_ind);
    ylabel(Outcome_Label)
    ylim(ylims_all); vline(0,'k')
    
    % Plotting peri-feedback period
    subplot( 1, 2, 2); hold on
    shadedErrorBar( ...
        psth_flanking_times(1:length(peri_feedback_mean_filt)),...
        peri_feedback_mean_filt,   ...
        peri_feedback_err_filt,    ...
        {trial_colors{trial_type_outcome_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_type_outcome_ind}(1)} );
    title({'Aligned to Feedback'})
    %ylabel(sprintf('Outcome Type %.0f',trial_group_ind))
    ylim(ylims_all); vline(0,'k')
    
    %% Targeted comparisons
    figure(comparison_full_fig)
    
    subplot(3,2,2*(trial_outcome_ind-1)+1); hold on
    plot(psth_flanking_times(1:length(peri_stim_mean_filt)),...
        peri_stim_mean_filt,  ...
        trial_colors{trial_type_outcome_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_type_outcome_ind})
    ylabel(trial_outcome_legends{trial_outcome_ind})
    vline(0)
    
    subplot(3,2,2*(trial_outcome_ind-1)+2); hold on
    plot( psth_flanking_times(1:length(peri_feedback_mean_filt)),peri_feedback_mean_filt,  trial_colors{trial_type_outcome_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_type_outcome_ind})
    ylabel(trial_outcome_legends{trial_outcome_ind})
    
    vline(0)
    
    
    
    %% Pairwise comparisons
    if do_all_pairwise_plots
        per
        subplot_row_inds = (trial_type_outcome_ind-1) * n_trial_type_vs_outcomes + [1:n_trial_type_vs_outcomes];
        subplot_col_inds = trial_type_outcome_ind + [0:n_trial_type_vs_outcomes:(n_trial_type_vs_outcomes * n_trial_type_vs_outcomes-1)];
        subplot_all_inds = unique([subplot_row_inds subplot_col_inds]);
        
        figure(pairwise_stim_fig)
        for subplot_ind = subplot_all_inds
            subplot( n_trial_type_vs_outcomes, n_trial_type_vs_outcomes , subplot_ind); hold on;
            %             shadedErrorBar( ...
            %                 psth_flanking_times, ...
            %                 peri_stim_mean,      ...
            %                 peri_stim_err,       ...
            %                 {trial_colors{trial_group_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_group_ind}} );
            
            plot(psth_flanking_times(1:length(peri_stim_mean_filt)),...
                peri_stim_mean_filt,  ...
                trial_colors{trial_type_outcome_ind},...
                'LineWidth',2,'markerfacecolor',trial_colors{trial_type_outcome_ind})
            Outcome_Label = sprintf('%s (%.0f)',trial_type_outcome_legend ,trial_type_outcome_ind);
            ylim(ylims_all); vline(0,'k')
            xlim(psth_flanking_xlim)
            
            if subplot_ind <= n_trial_type_vs_outcomes
                title(sprintf('%s (%.0f)', Outcome_Label, peri_feedback_n))
                set(gca,'YTick',[])
                set(gca,'XTick',[])
            elseif mod(subplot_ind,2) == 1 && (floor(subplot_ind/n_trial_type_vs_outcomes) == mod(subplot_ind,n_trial_type_vs_outcomes))
                axis tight
            else
                axis off
            end
        end
        
        figure(pairwise_feedback_fig)
        for subplot_ind = subplot_all_inds
            
            % Plotting peri-feedback period
            subplot( n_trial_type_vs_outcomes, n_trial_type_vs_outcomes , subplot_ind); hold on;
            %             shadedErrorBar( ...
            %                 psth_flanking_times,...
            %                 peri_feedback_mean,   ...
            %                 peri_feedback_err,    ...
            %                 {trial_colors{trial_group_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_group_ind}} );
            plot(psth_flanking_times(1:length(peri_feedback_mean)),...
                peri_feedback_mean,  ...
                trial_colors{trial_type_outcome_ind},...
                'LineWidth',2,'markerfacecolor',trial_colors{trial_type_outcome_ind})
            
            %ylabel(sprintf('Outcome Type %.0f',trial_group_ind))
            ylim(ylims_all); vline(0,'k')
            xlim(psth_flanking_xlim)
            
            if subplot_ind <= n_trial_type_vs_outcomes
                title(sprintf('%s (%.0f)', Outcome_Label, peri_feedback_n))
                set(gca,'YTick',[])
                set(gca,'XTick',[])
            elseif  mod(subplot_ind,2) == 1 && ...
                    (floor(subplot_ind/n_trial_type_vs_outcomes) == mod(subplot_ind,n_trial_type_vs_outcomes))
                axis tight
                % Leave these axes alone.
            else
                % Turn off these axes for readability.
                axis off
            end
            
        end
    end
    
end

end
tilefigs
disp('Done')
beep;pause(0.5); beep










%%



%%

