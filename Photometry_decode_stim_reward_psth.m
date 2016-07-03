
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

TODO ::
- Probably worth using FilterX compiled version
- filters can be combined in transfer function form.


%}

%% MISC DEBUGGING CODE




%%
dbstop if error

%file_directory = '/Volumes/Macintosh HD/Users/avaughan/Dropbox/KepecsLab (1)/_Alex/Photometry/bpod data/';
%file_directory = '~/Documents/MATLAB/_projects_CSHL/photometry_analysis/bpod_data/';
file_directory = '/Users/avaughan/Dropbox/MATLAB/_projects_CSHL/Photometry/data/';

file_names = {};
trials_to_skip = [];
trial_outcomes = {};
notch_frequency = NaN;


% Sep 01
% At some point before this, audio was unplugged :|
% However, we can see audio in channel 3 for confirmation.

data_to_analyze = 'Lockin_test'

switch data_to_analyze
    
     case 'GCaMP1'
        % GCaMP1 (black) - needs
        % WHERE IS SEP 01?
        file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Sep02_2015_Session1.mat';
        %file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Sep03_2015_Session1.mat';  % Very noisy
        file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Sep10_2015_Session1.mat';  % --> Very weird.
        file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Sep11_2015_Session1.mat';
        notch_frequency = 13; % Some 13hz for Sep02/03, but stronger 20hz for Sep10.
        
        %file_names{end+1} = 'VIP_Ach_black_SO_Training_NIDAQ_Sep01_2015_Session1.mat';     % AUDIO OK
    
    case 'GCaMP2'
        % GCaMP2 (pink) - needs 20hz filter
        file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Sep01_2015_Session2.mat'; % Good session #1
        file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Sep02_2015_Session1.mat'; % Good session #2
        % file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Sep03_2015_Session1.mat'; % very bad noise
        % file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Sep10_2015_Session2.mat'; % very bad noise
        notch_frequency = 20;
        
    case 'GCaMP3'
        % Lockin signal was wrong for Sep 24th and most of Sep 25.
        % Should be good for Sep 28+ %Also, lockin freq. changed to 1000hz.
        %file_names{end+1} = 'VIP_GCaMP_3_pink_SO_Training_NIDAQ_Sep23_2015_Session3.mat'; % Bad lockin
        file_names{end+1} = 'VIP_GCaMP_3_pink_SO_Training_NIDAQ_Sep24_2015_Session2.mat'; % Maybe punishment?
        file_names{end+1} = 'VIP_GCaMP_3_pink_SO_Training_NIDAQ_Sep29_2015_Session2.mat'; % Maybe punishment
        file_names{end+1} = 'VIP_GCaMP_3_pink_SO_Training_NIDAQ_Sep30_2015_Session2.mat'; % Maybe punishment?
        %notch_frequency = 20;
        
        
    case 'GCaMP4'
        
        % Generally useless.
        %file_names{end+1} = 'VIP_GCaMP_4_black_SO_Training_NIDAQ_Sep22_2015_Session1.mat' % Crap        Crap
        %file_names{end+1} = 'VIP_GCaMP_4_black_SO_Training_NIDAQ_Sep22_2015_Session2.mat' % Crap
        %file_names{end+1} = 'VIP_GCaMP_4_black_SO_Training_NIDAQ_Sep23_2015_Session1.mat' % Crap
        
    case 'GCaMP5'
        %
        %file_names{end+1} = 'VIP_GCaMP_5_black_SO_Training_NIDAQ_Oct02_2015_Session1.mat'
        %file_names{end+1} = 'VIP_GCaMP_5_black_SO_Training_NIDAQ_Oct03_2015_Session1.mat' % --> low trial #
        %file_names{end+1} = 'VIP_GCaMP_5_black_SO_Training_NIDAQ_Oct04_2015_Session1.mat' % -->
        %file_names{end+1} = 'VIP_GCaMP_5_black_SO_Training_NIDAQ_Oct05_2015_Session1.mat'         Crappy response
        file_names{end+1} = 'VIP_GCaMP5Co5_SO_Training_NIDAQ_Oct08_2015_Session1.mat'
        %file_names{end+1} = 'VIP_GCaMP_5_black_SO_Training_NIDAQ_Oct09_2015_Session1.mat'  % Bad anticipatory licking.
        %file_names{end+1} = VIP_GCaMP_5_black_SO_Training_NIDAQ_Oct13_2015_Session3.mat'   % Crappy response
        
    case 'Chat_pink'
        
        % Chat files.
        file_names{end+1} = 'Chat_GCaMP_pink_SO_Training_NIDAQ_Sep25_2015_Session2.mat';
        file_names{end+1} = 'Chat_GCaMP_pink_SO_Training_NIDAQ_Sep28_2015_Session2.mat';
        file_names{end+1} = 'Chat_GCaMP_pink_SO_Training_NIDAQ_Sep29_2015_Session1.mat';
        
        
    case 'GCaMP1 + GCaMP2 + GCaMP3'
        
        file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Sep02_2015_Session1.mat';
        file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Sep10_2015_Session1.mat';  % --> Very weird.
        file_names{end+1} = 'VIP_GCaMP1_B124_SO_Training_NIDAQ_Sep11_2015_Session1.mat';
        file_names{end+1} = 'VIP_GCaMP_3_pink_SO_Training_NIDAQ_Sep24_2015_Session2.mat'; % Maybe punishment?
        file_names{end+1} = 'VIP_GCaMP_3_pink_SO_Training_NIDAQ_Sep29_2015_Session2.mat'; % Maybe punishment
        file_names{end+1} = 'VIP_GCaMP_3_pink_SO_Training_NIDAQ_Sep30_2015_Session2.mat'; % Maybe punishment?
        file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Sep01_2015_Session2.mat'; % Good session #1
        file_names{end+1} = 'VIP_GCaMP_2_SO_Training_NIDAQ_Sep02_2015_Session1.mat'; % Good session #2

    case 'Lockin_test'
        %file_names{end+1} = 'VIP_GCaMP_3_pink_SO_Training_NIDAQ_Sep30_2015_Session2.mat'; % Maybe punishment?
        file_names{end+1} = 'VIP_GCaMP5Co5_SO_Training_NIDAQ_Oct08_2015_Session1.mat'

    otherwise
        % Lockin 
        error('Choose some files, bro.')
        
end
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

trial_outcomes        = {  [0 4],       [1   2],    [3 5]       };
trial_outcome_fields  = { 'Punish',     'Reward', 'PostUS'    };
trial_outcome_legends = { 'Punishment', 'Reward', 'Omission'  };
trial_colors          = {  'm','r' ,    'g','c',  'b','k'           };
trial_type_plot_order = [1 2 4 3 5 6]; % For convenience with what's on top;

%%%trial_colors          = { 'r',          'g',      'b'           };

% trial_outcomes = { 0,4,1,2,3,5 };
% trial_feedback_fields = { 'Punish', 'Punish',    'Reward','Reward', 'PostUS', 'PostUS'      };
% trial_legends         = { 'Punish (pre-lick)', 'Punish (no-lick)',    'Reward (pre-lick)','Reward (no-lick)', 'Omission (pre-lick)', 'Omission (no-lick)'      };
% trial_colors          = { 'm', 'r' ,  'c','g',      'b','k'           };

% Trial types are 1 (predominantly reward) and 0 (predominantly punishment)
trail_stimuli        = [1 2]; % <-- drawn from data below, but this is usually what
trial_type_legends = {'E(Reward)','E(Punish)'};
%it is.


%% Analysis parameters

% DO ANALYSIS OR JUST PRETEND
actually_do_analysis            = 1;

% COMBINE ALL FILES ON A TRIAL-BY-TRIAL BASIS?
combine_trials_across_sessions  = 1;
%normalize_each_session          = 1;

% Do all plots, or only sanity-check plots?
do_all_plots                    = 1;

% Bleaching and trial-length figures
do_sanity_figs  = 0;

% Do all pairwise plots?  This can get annoying.
do_all_pairwise_plots           = 0;

% Downsample GCaMP trials?  10 --> downsample 10x.
downsample_int                  = 1;

% Lockin frequency (hz)
lockin_frequency = 470;

% Zeroing function (as a function of sampling rate, which can vary based on
% downsampling.
zeroing_fn = @(sample_rate) sample_rate/2;

% Exclude trials with high variance?
% Excludes trials where var(trial) > [ variance_threshold * std(var(all_trials))
variance_threshold = 3;

% SKIP SOME TRIALS FOR ALL FILES?
% Only Trials 5:55
%trials_to_skip_for_all_fn = @(n_trials) [1:5 56:n_trials n_trials-[0:4]];
% Skip 1st 5 and last 5 trials, and anything over trial 250
trials_to_skip_for_all_fn = @(n_trials) [1:5, n_trials-[0:4], 251:1000];

% Filtering?
do_filtering = 0;
filter_type = {'notch'}; % 'notch' or 'lowpass', 'bandpass'

% Time (s) on the either side of the alignment cue
% _range is the full range for filtering, etc.
% _display is the actual xlims on display
psth_flanking_range   = [-2.0 2.1]; % max pre is ~2s because of stimulus presentation.  Can be changed though
psth_flanking_xlim    = [-1.5 2];

% Y range for all figures - defined in units of STDEV of all samples of the data
ylims_GCaMP_fn = @(all_GCaMP) nanstd(all_GCaMP(:)) * [-2 2];
ylims_audio_fn = @(all_audio) nanstd(all_audio(:)) * [-1.5 1.5];

% Plot parameters
patch_alpha = 0.9;

%%
assert(~isempty(file_names),'file_names is empty :: You should probably select some files to analyze')


sample_rate = 10000 / downsample_int;
% Size of window for zeroing response.
zeroing_samples = zeroing_fn(sample_rate);

% FILTERING
% Low pass
Fpass = 45; Fstop = 50; Ap = 2; Ast = 20;
lowpass.d = designfilt(          ...
    'lowpassiir',                ...
    'PassbandFrequency',  Fpass, ...
    'StopbandFrequency',  Fstop, ...
    'PassbandRipple',     Ap,    ...
    'StopbandAttenuation',Ast,   ...
    'SampleRate',         sample_rate );
[lowpass.b,lowpass.a] = tf(lowpass.d);

% Notch (bandstop)
if isnan(notch_frequency)
    notch_frequency = 20; % or 34...
end
notch.d = designfilt('bandstopiir',...
    'FilterOrder',10,              ...
    'HalfPowerFrequency1',notch_frequency * 0.9, ...
    'HalfPowerFrequency2',notch_frequency * 1.1, ...
    'SampleRate',sample_rate );
[notch.b,notch.a] = tf(notch.d);

% Bandpass.
bandpass.d = designfilt(...
    'bandpassiir',...
    'FilterOrder',         10,     ...
    'HalfPowerFrequency1', 1,    ...
    'HalfPowerFrequency2', 30,     ...
    'SampleRate',          sample_rate );
[bandpass.b,bandpass.a] = tf(bandpass.d);

%%
switch do_filtering
    case 1
        filter_string = ' :: FILTERING ON';
    case 0
        filter_string = ' :: FILTERING OFF';
end


for i_file = 1:length(file_names);
    
    % Load file
    file_name = file_names{i_file};
    load([file_directory file_name]);
    
    fprintf('\nLoaded file %s\n',file_name)
    
    %% Trial params
    n_trials = length(SessionData.TrialTypes);
    n_channels = size(SessionData.NidaqData{1},2);
    fluorescent_channel= 1; % For Ai0
    raw_camera_channel = 2;
    audio_channel      = 3;
    
    %% PSTH samples.
    psth_flanking_times   = linspace( psth_flanking_range(1),psth_flanking_range(2), sample_rate*diff(psth_flanking_range)+1  );
    psth_flanking_samples = psth_flanking_times * sample_rate;
    
    %% Pre-process NidaqData by downsampling fluorescent channel
    for i_trial = 1:n_trials
        %SessionData.NidaqData{i_trial} = medfilt1(SessionData.NidaqData{i_trial}(1:downsample_int:end,:),ceil(downsample_int/2),'truncate');
        SessionData.NidaqData{i_trial} = SessionData.NidaqData{i_trial}(1:downsample_int:end,:);
    end
    
    %% Pre-process NidaqData by tossing trials with absurdly high variance
    all_trial_variances = zeros(length(SessionData.NidaqData),1);
    for i_trial = 1:n_trials
        all_trial_variances(i_trial) = var(SessionData.NidaqData{i_trial}(:,fluorescent_channel));
    end
    trials_with_excess_variance = find(abs(all_trial_variances) > (variance_threshold*std(all_trial_variances)));
    
    fprintf('  Excluding %.0f trials (%.0f%%) b/c of variance > %.1f sigma.\n',length(trials_with_excess_variance),100*length(trials_with_excess_variance)/length(SessionData.NidaqData),variance_threshold)
  
    
  
    
    %% Default if not specified above.
    if isempty(trial_outcomes)
        trial_outcomes = num2cell(unique(SessionData.TrialOutcome));
    end
    
    %% Generate nd_trial_typees and nd_trial_outcomes to iterate over.
    trail_stimuli    = unique(SessionData.TrialTypes);
    
    % Index
    [nd_trail_stimuli,nd_outcome_types] = ndgrid(1:length(trail_stimuli),1:length(trial_outcomes));
    trial_type_matrix = [nd_trail_stimuli(:) nd_outcome_types(:)];
    n_trial_types = size(trial_type_matrix,1);
    
    % Trial sizes etc.
    [trial_sizes, postUS_times] = deal([]);
    for i = 1:SessionData.nTrials,
        trial_sizes(i) = size(SessionData.NidaqData{i},1)/sample_rate;
    end
    for i = 1:SessionData.nTrials,
        try,
            % SO (stimulus-outcome_ training
            postUS_times(i) = SessionData.RawEvents.Trial{i}.States.PostUS(1);
        catch
            % AudGoNoGo task.
            postUS_times(i) = SessionData.RawEvents.Trial{i}.States.WaitForLick(2);
        end
    end
    trial_length_cutoff = find(postUS_times > trial_sizes,1);

    
    %%
    fprintf('Available trial types: %s\n', mat2str(cell2mat(trial_outcomes)))
    fprintf('Available trial outcomes: %s\n',mat2str(cell2mat(trial_outcomes)))
    fprintf(' %.0f pairwise combinations\n',n_trial_types)
    
    if do_sanity_figs
        
        %% Trial Size Figure and cutoff for short trials
        namedFigure(sprintf('Trial Sizes :: %s',file_name)); clf; hold on
        plot(trial_sizes,'b')
        plot(postUS_times,'g')
        plot(postUS_times > trial_sizes,'r.')
        legend({'Recording Length','PostUS Time','MISMATCH'})
        
        if isempty(trial_length_cutoff), trial_length_cutoff = SessionData.nTrials+1; end
        fprintf('trial_length_cutoff is set to trial #%.0f\n',trial_length_cutoff)
        
        
        
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
        
    end
    
    
    %%  Loop across trial types and aggregate data
    sample_trials = cell(size(n_trial_types));
    for i_trial_type  = 1:size(trial_type_matrix,1)
        
        %% Aggregate trial data
        % Indices corresponding to appropriate stimulus / outcome.
        i_trial_stimulus            = trial_type_matrix(i_trial_type,1);
        i_trial_outcome             = trial_type_matrix(i_trial_type,2);
        this_trial_stimulus         = trail_stimuli(i_trial_stimulus);
        this_trial_outcome          = trial_outcomes{i_trial_outcome};
        this_trial_outcome_field    = trial_outcome_fields{i_trial_outcome}; % 'Punish' or 'DeliverPunish'
        this_trial_type_legend      = sprintf('%s > %s',trial_type_legends{i_trial_stimulus},trial_outcome_legends{i_trial_outcome});
        
        fprintf('%s\n',this_trial_type_legend)
        fprintf('  trial_type %.0f, trial_outcome %s ::',this_trial_stimulus,mat2str(this_trial_outcome))
        
        % Collect trials (either grouping more than one
        if length(this_trial_outcome) > 1
            % Group more than one condition into the same bin
            trials_this_trial_type = find(sum( bsxfun(@eq,SessionData.TrialOutcome,this_trial_outcome')));
        else
            % Find only a single condition
            trials_this_trial_type = find(bsxfun(@eq,SessionData.TrialOutcome,this_trial_outcome'));
        end
        fprintf(' %03.0f trials\n',length(trials_this_trial_type))
        % Cull by stimulus type
        trials_this_trial_type = intersect(trials_this_trial_type,find(SessionData.TrialTypes == this_trial_stimulus));
        fprintf('  culling by trial_stimulus = %.0f\t    :: %03.0f trials\n',this_trial_stimulus,length(trials_this_trial_type))
        % Cull by trial_length_cutoff
        trials_this_trial_type  = trials_this_trial_type(trials_this_trial_type < trial_length_cutoff);
        fprintf('  culling by trial_length_cutoff :: %03.0f trials\n',length(trials_this_trial_type))
        % Cull by trials_to_skip
        trials_this_trial_type = setdiff(trials_this_trial_type,trials_to_skip);
        trials_this_trial_type = setdiff(trials_this_trial_type,trials_to_skip_for_all_fn(n_trials));
        trials_this_trial_type = round(trials_this_trial_type);
        fprintf('  culling by trials_to_skip\t    :: %03.0f trials.\n',length(trials_this_trial_type))
        trials_this_trial_type = setdiff(trials_this_trial_type,trials_with_excess_variance);
        fprintf('  culling by high variance \t    :: %03.0f trials.\n',length(trials_this_trial_type))
        
        % Do we have enough trials?
        n_this_trial_type = length(trials_this_trial_type);
        if n_this_trial_type < 2
            fprintf('*** WARNING :: <1 TRIAL -- SKIPPING ANALYSIS ***\n')
            continue
        end
        
        % Skip analysis for debugging?
        if ~actually_do_analysis
            continue
        end
        
        % Initialize holding variables.
        [ peri_stimulus_GCaMP{i_file,i_trial_type}, ...
            peri_feedback_GCaMP{i_file,i_trial_type}, ...
            peri_stimulus_GCaMP_processed{i_file,i_trial_type}, ...
            peri_feedback_GCaMP_processed{i_file,i_trial_type},...
            peri_stimulus_audio{i_file,i_trial_type},...
            peri_feedback_audio{i_file,i_trial_type} ] =  deal(nan(n_this_trial_type,length(psth_flanking_samples)));
        [ peri_stimulus_licks{i_file,i_trial_type}, ...
            peri_feedback_licks{i_file,i_trial_type} ] = deal(sparse(n_this_trial_type,length(psth_flanking_samples)));
        
        %% Loop over each trial and extract signal
        for i_trial = 1:length(trials_this_trial_type),
            trial = trials_this_trial_type(i_trial);
            
            %% Pre-process NidaqData by decoding
            
            do_lockin_plot = 0;
            tic
            keyboard
            [smoothed_output, demodulation_freq, demodulation_phase] = decode_lockin_fn( ...
                SessionData.NidaqData{trial}( : , raw_camera_channel),...  % To decode
                [] , ... % Lockin reference sinusoid
                SessionData.NidaqData{trial}( : , fluorescent_channel),... % Reference.
                lockin_frequency, sample_rate, ...
                do_lockin_plot); % do plot?
            toc
            %  [smoothed_output, demodulation_freq, demodulation_phase] = decode_lockin( ...
            %  SessionData.NidaqData{trial}( stimulus_t_indices(1:trial_end_sample) , raw_camera_channel),...  % To decode
            %  [] , ... % Lockin reference sinusoid
            %  SessionData.NidaqData{trial}( stimulus_t_indices(1:trial_end_sample) , fluorescent_channel),... % Reference.
            %  lockin_frequency, sample_rate, ...
            %  do_lockin_plot); % do plot?
            
            figure(1); %clf
            % subplot(2,1,1)
            % pwelch(SessionData.NidaqData{trial}(:,raw_camera_channel), 1024,512,[],sample_rate)
            subplot(2,2,i_trial); cla; hold on
            h = plot(zscore(smoothed_output(sample_rate/2:end-sample_rate/2) ));
            h(2) = plot(zscore(SessionData.NidaqData{trial}( sample_rate/2:end-sample_rate/2,fluorescent_channel)) );
            set(h,'LineWidth',2)
            if i_trial == 1
                legend('Hardware','Software')
            end
            axis tight
            
            if i_trial == 4
                beep; pause(0.5); beep
                keyboard
            end
            continue
            
            % Replace in-place to make things simple downstream
            SessionData.NidaqData{trial}( : , raw_camera_channel) = smoothed_output;
            
            %% Peri-Stimulus period
            % GCaMP response
            stimulus_start_time = SessionData.RawEvents.Trial{trial}.States.DeliverStimulus(1);
            % Time index for stimulus display period
            stimulus_t_indices  = round( stimulus_start_time * sample_rate + psth_flanking_samples(:) );
            available_samples   = size( SessionData.NidaqData{trial} , 1) - stimulus_start_time * sample_rate;
            if available_samples < 0,
                % Leave this row as NaNs.
                fprintf('available_samples < 0 for some reason (%.0f) :: trial %.0f :: stimulus\n',available_samples,trial)
                continue
            end
            trial_end_sample = round(min(length(psth_flanking_samples),available_samples));
            
            % Symmetric window around start
            % stimulus_zeroing_window = (stimulus_start_time * sample_rate) + round(linspace(-floor(zeroing_samples/2),floor(zeroing_samples/2)-1,zeroing_samples));
            % Causal window around start
            stimulus_zeroing_window = (stimulus_start_time * sample_rate) + round(linspace(-zeroing_samples,-1,zeroing_samples));
            
            % Populate peri_stimulus_GCaMP matrix with response normalized by mean around t=0;
            peri_stimulus_GCaMP{i_file,i_trial_type}( i_trial, 1:trial_end_sample ) = ...
                SessionData.NidaqData{trial}( stimulus_t_indices(1:trial_end_sample) , fluorescent_channel)' ...
                - mean(SessionData.NidaqData{trial}( round(stimulus_zeroing_window) , fluorescent_channel));
            peri_stimulus_audio{i_file,i_trial_type}( i_trial, 1:trial_end_sample ) = ...
                mean(SessionData.NidaqData{trial}( stimulus_t_indices(1:trial_end_sample) , audio_channel ),2)' ...
                - mean(mean(SessionData.NidaqData{trial}( round(stimulus_zeroing_window) , audio_channel),2));
           
            
            
            %% Peri-Feedback period
            
            % Find time for beginning of feedback.
            feedback_time       = SessionData.RawEvents.Trial{trial}.States.(this_trial_outcome_field);
            feedback_start_time = feedback_time(1);
            assert( ~isnan(feedback_start_time) , 'Feedback start time is a NaN - usually a sign that you''re looking at the wrong behavioral timepoint.')
            
            % Time index for feedback samples.
            feedback_t_indices = round( feedback_start_time * sample_rate + psth_flanking_samples(:) );
            available_samples = size( SessionData.NidaqData{trial} , 1) - feedback_start_time * sample_rate;
            if available_samples < 0,
                fprintf('available_samples < 0 for some reason :: trial %.0f :: reward\n',trial)
                continue
            end
            trial_end_sample = round(min(length(psth_flanking_samples),available_samples));
            
            % Normalizing by feedback start time to be dF/F is not very
            % helpful the baseline value is below zero.  So we just do dF.
            feedback_zeroing_window = round(feedback_start_time * sample_rate + linspace(-zeroing_samples,-1,zeroing_samples));
            
            peri_feedback_GCaMP{i_file,i_trial_type}( i_trial, 1:trial_end_sample ) = ...
                ( SessionData.NidaqData{trial}( feedback_t_indices(1:trial_end_sample) , fluorescent_channel)' ...
                - mean(SessionData.NidaqData{trial}( feedback_zeroing_window , fluorescent_channel)));
            
            peri_feedback_audio{i_file,i_trial_type}( i_trial, 1:trial_end_sample ) = ...
                mean( SessionData.NidaqData{trial}( feedback_t_indices(1:trial_end_sample) , audio_channel ),2)' ...
                - mean(mean(SessionData.NidaqData{trial}( feedback_zeroing_window , audio_channel),2));
            %/ SessionData.NidaqData{trial}( round(feedback_start_time * sample_rate) , fluorescent_channel) ;
            
            %% Filtering
            
            % Notch or smoothing filtering.
            % Make sure we only filter the parts of this matrix that are
            % sensible - outside of 1:end_sample may be nan at this point.
            
            if do_filtering
                % Do Filtering for both Stimulus and Feedback signals - via
                % dummy variables for brevity.
                input_signals  = { peri_stimulus_GCaMP{i_file,i_trial_type} ,           peri_feedback_GCaMP{i_file,i_trial_type}           };
                output_signals = { peri_stimulus_GCaMP_processed{i_file,i_trial_type} , peri_feedback_GCaMP_processed{i_file,i_trial_type} };
                
                for i_signal = [1 2]
                    output_signals{i_signal}(i_trial,1:trial_end_sample) = input_signals{i_signal}(end,1:trial_end_sample);
                    if any(strcmp(filter_type,'bandpass'))
                        output_signals{i_signal}( end, 1:trial_end_sample ) = filtfilt(bandpass.d, output_signals{i_signal}(end, 1:trial_end_sample)');
                        %peri_stimulus_GCaMP_processed{trial_type_outcome_ind}( end, 1:end_sample ) = FiltFiltM(bandpass.b, bandpass.a, peri_stimulus_GCaMP_processed{trial_type_outcome_ind}(end, 1:end_sample)');
                    end
                    if any(strcmp(filter_type,'lowpass'))
                        output_signals{i_signal}( end, 1:trial_end_sample )  = filtfilt(lowpass.d, output_signals{i_signal}(end, 1:trial_end_sample)');
                        %peri_stimulus_GCaMP_processed{trial_type_outcome_ind}( end, 1:end_sample ) = FiltFiltM(lowpass.b, lowpass.a, peri_stimulus_GCaMP_processed{trial_type_outcome_ind}(end, 1:end_sample)');
                    end
                    if any(strcmp(filter_type,'notch'))
                        output_signals{i_signal}( end, 1:trial_end_sample )  = filtfilt(notch.d, output_signals{i_signal}(end, 1:trial_end_sample)');
                        %peri_stimulus_GCaMP_processed{trial_type_outcome_ind}( end, 1:end_sample ) = FiltFiltM(notch.b, notch.a, peri_stimulus_GCaMP_processed{trial_type_outcome_ind}(end, 1:end_sample)');
                    end
                end
                peri_stimulus_GCaMP_processed{i_file,i_trial_type} = output_signals{1};
                peri_feedback_GCaMP_processed{i_file,i_trial_type} = output_signals{2};
            else
                % No filtering.
                peri_stimulus_GCaMP_processed{i_file,i_trial_type} = peri_stimulus_GCaMP{i_file,i_trial_type};
                peri_feedback_GCaMP_processed{i_file,i_trial_type} = peri_feedback_GCaMP{i_file,i_trial_type};
            end
            
            
            %% Calculate lick rate.
            % Stimulus
            stimulus_start_time = SessionData.RawEvents.Trial{trial}.States.DeliverStimulus(1);
            available_samples = size( SessionData.NidaqData{trial} , 1) - stimulus_start_time * sample_rate;
            
            % Gather licks
            if isfield(SessionData.RawEvents.Trial{trial}.Events,'Port2In')
                all_licks = SessionData.RawEvents.Trial{trial}.Events.Port2In;
                %fprintf('    Trial %.0f :: %.0f Port2In events\n',trial,length(all_licks ))
            else
                %fprintf('    Trial %.0f :: No Port2In events!\n',trial)
                all_licks = [];
            end
            
            % Find licks in stimulus/feedback window.
            stimulus_window = stimulus_start_time + psth_flanking_range;
            feedback_window = feedback_start_time + psth_flanking_range;
            % Find licks and align to beginning of stimulus/feedback window
            this_trial_peri_stimulus_licks = all_licks( stimulus_window(1) <= all_licks & all_licks <= stimulus_window(2)) - stimulus_window(1);
            this_trial_peri_feedback_licks = all_licks( feedback_window(1) <= all_licks & all_licks <= feedback_window(2)) - feedback_window(1);
            % Convert to samples after beginning of stimulus/feedback window
            peri_stimulus_licks{i_file,i_trial_type}(i_trial, ceil( sample_rate*( this_trial_peri_stimulus_licks ))) = 1;
            peri_feedback_licks{i_file,i_trial_type}(i_trial, ceil( sample_rate*( this_trial_peri_feedback_licks ))) = 1;
            
        end % END Iterating over trials.
        
        %% Error handling on the length of the peri_feedback_GCaMP etc.
        if size(peri_feedback_GCaMP_processed{i_file,i_trial_type},2)  < length(psth_flanking_samples)
            fprintf('  *** WARNING *** :: trial_type_outcome_ind %.0f :: expected %.0f samples but only found %.0f\np',...
                i_trial_type,...
                length(psth_flanking_samples),...
                size(peri_feedback_GCaMP_processed{i_file,i_trial_type},2))
        end
        
    end % END TRIAL TYPES (ANALYSIS)
    
end % Files (ANALYSIS)

if ~actually_do_analysis
    return
end

%% Prepare for plotting

% Index for looping through the plotting code.  If combining all
% sessions, then pass a matrix of file indices.  Otherwise, just pass one
% the index for a single fie.
if combine_trials_across_sessions
    % Passed as matrix - ie., plot all at once.
    files_to_plot = {1:length(file_names)};
else
    % Passed as single values - ie., plot separately.
    files_to_plot = num2cell(1:length(file_names));
end

% Generate consistent y limits across all sessions / trial types.
all_GCaMP = vertcat(peri_stimulus_GCaMP{:},peri_feedback_GCaMP{:});
all_audio = vertcat(peri_stimulus_audio{:},peri_feedback_audio{:});
ylims_GCaMP = ylims_GCaMP_fn(all_GCaMP);
ylims_audio = ylims_audio_fn(all_audio);
clear all_GCaMP all_audio;


%% PLOTTING
for plotting_file_inds = files_to_plot
    plotting_file_inds = plotting_file_inds{1};
    
    %% If skipping plotting, continue looping here.
    if combine_trials_across_sessions
        if ~(i_file == length(file_names))
            fprintf('Skipping plots for now - combining later.\n')
            continue
        else
            figure_prefix = data_to_analyze;
        end
    else
        figure_prefix = file_names{plotting_file_inds};
    end
    
    % Initialize figures.
    audio_fig                   = namedFigure([figure_prefix '__audio' filter_string]); clf
    if do_all_plots
        average_fig                 = namedFigure([figure_prefix '__average' filter_string]); clf
        all_trials_fig              = namedFigure([figure_prefix '__all_trials' filter_string]); clf
        spectrogram_fig             = namedFigure([figure_prefix '__spectrogram' filter_string]); clf
        combined_fig                = namedFigure([figure_prefix '__combined' filter_string]); clf
        sample_fig                  = namedFigure([figure_prefix '__sample' filter_string]); clf
        pairwise_stimulus_fig           = namedFigure([figure_prefix '__pairwise_stim' filter_string]); clf
        pairwise_feedback_fig       = namedFigure([figure_prefix '__pairwise_feedback' filter_string]); clf
        comparison_full_fig         = namedFigure([figure_prefix '__comparison_figure' filter_string]); clf
        lick_raster_fig             = namedFigure([figure_prefix '__lick_raster' filter_string]); clf
        if do_all_pairwise_plots
            pairwise_stimulus_fig       = namedFigure([figure_prefix '__pairwise_stim' filter_string]); clf
            pairwise_feedback_fig   = namedFigure([figure_prefix '__pairwise_feedback' filter_string]); clf
        end
    end
    
    
    %% PLOTTING
    for i_trial_type  = trial_type_plot_order
        
        %%
        i_trial_stimulus = trial_type_matrix(i_trial_type,1);
        i_trial_outcome  = trial_type_matrix(i_trial_type,2);
        this_trial_stimulus     = trail_stimuli(i_trial_stimulus);
        this_trial_outcome      = trial_outcomes{i_trial_outcome};
        this_trial_type_legend = sprintf('%s > %s',trial_type_legends{i_trial_stimulus},trial_outcome_legends{i_trial_outcome});
        trail_type_outcome_yaxis = {trial_type_legends{i_trial_stimulus},trial_outcome_legends{i_trial_outcome}};
        fprintf('PLOTTING :: %s\n',this_trial_type_legend)
        
        %% Prep for plotting
        
        % Cell - subset of the full cell containing sesions x trail_stimuli
        this_peri_stimulus_GCaMP           = {peri_stimulus_GCaMP{plotting_file_inds,i_trial_type}};
        this_peri_feedback_GCaMP           = {peri_feedback_GCaMP{plotting_file_inds,i_trial_type}};
        this_peri_stimulus_GCaMP_processed = {peri_stimulus_GCaMP_processed{plotting_file_inds,i_trial_type}};
        this_peri_feedback_GCaMP_processed = {peri_feedback_GCaMP_processed{plotting_file_inds,i_trial_type}};
        this_peri_stimulus_audio           = {peri_stimulus_audio{plotting_file_inds,i_trial_type}};
        this_peri_feedback_audio           = {peri_feedback_audio{plotting_file_inds,i_trial_type}};
        this_peri_stimulus_licks           = {peri_stimulus_licks{plotting_file_inds,i_trial_type}};
        this_peri_feedback_licks           = {peri_feedback_licks{plotting_file_inds,i_trial_type}};

        % Matrix - all trials across sessions x file_types
        this_peri_stimulus_GCaMP_all_trials           = vertcat(this_peri_stimulus_GCaMP{:});
        this_peri_feedback_GCaMP_all_trials           = vertcat(this_peri_feedback_GCaMP{:});
        this_peri_stimulus_GCaMP_processed_all_trials = vertcat(this_peri_stimulus_GCaMP{:});
        this_peri_feedback_GCaMP_processed_all_trials = vertcat(this_peri_feedback_GCaMP{:});
        this_peri_stimulus_audio_all_trials           = vertcat(this_peri_stimulus_audio{:});
        this_peri_feedback_audio_all_trials           = vertcat(this_peri_feedback_audio{:});
        this_peri_stimulus_licks_all_trials           = vertcat(this_peri_stimulus_licks{:});
        this_peri_feedback_licks_all_trials           = vertcat(this_peri_feedback_licks{:});
        
        n_this_trial_type = size(this_peri_stimulus_GCaMP_all_trials,1);

        % Why?
        %this_peri_stimulus_GCaMP_all_trials(            isnan(this_peri_stimulus_GCaMP_all_trials))           = 0;
        %this_peri_stimulus_GCaMP_processed_all_trials(  isnan(this_peri_stimulus_GCaMP_processed_all_trials)) = 0;
        %this_peri_feedback_GCaMP_all_trials(            isnan(this_peri_feedback_GCaMP_all_trials))           = 0;
        %this_peri_feedback_GCaMP_processed_all_trials(  isnan(this_peri_feedback_GCaMP_processed_all_trials)) = 0;
        
        % This two-step averaging seems a little complicated, but it allows for cross-session
        % averaging in two steps .  Note that we use the "processed" version
        % for plotting by default - if filtering is off, these are
        % standins.
        % Average standard error = sqrt( sum( err^2 + err^2 + ...)) - assumes no covariance, which is somewhat cheating.
        
        % STIM :: Mean and error within each session
        err_fn = @(x) nanstd(x) ./ sqrt(n_this_trial_type);
        peri_stimulus_GCaMP_processed_session_means = cellfun(@nanmean,this_peri_stimulus_GCaMP_processed,'UniformOutput',false);
        peri_feedback_GCaMP_processed_session_means = cellfun(@nanmean,this_peri_feedback_GCaMP_processed,'UniformOutput',false);
        peri_stimulus_GCaMP_processed_session_errs  = cellfun(err_fn , this_peri_stimulus_GCaMP_processed,'UniformOutput',false);
        peri_feedback_GCaMP_processed_session_errs  = cellfun(err_fn , this_peri_feedback_GCaMP_processed,'UniformOutput',false);
        peri_stimulus_GCaMP_processed_mean = nanmean( vertcat(peri_stimulus_GCaMP_processed_session_means{:}) , 1 );
        peri_feedback_GCaMP_processed_mean = nanmean( vertcat(peri_feedback_GCaMP_processed_session_means{:}) , 1 );
        peri_stimulus_GCaMP_processed_err  = sqrt(nansum(vertcat(peri_stimulus_GCaMP_processed_session_errs{:}).^2,1));
        peri_feedback_GCaMP_processed_err  = sqrt(nansum(vertcat(peri_feedback_GCaMP_processed_session_errs{:}).^2, 1));
        
        % FEEDBACK :: Mean/standard error within each session
        
        % Find representative sample (ie., most like the average);
        %peri_feedback_GCaMP_processed_nanafe = peri_feedback_GCaMP_processed{trial_group_ind};
        %peri_feedback_GCaMP_processed_nanafe(isnan(peri_feedback_GCaMP_processed_nanafe )) = 0;
        template_sum_squared_error = sum(bsxfun(@minus,this_peri_feedback_GCaMP_all_trials,nanmean(this_peri_feedback_GCaMP_all_trials,1)).^2,2);
        this_sample_trial = find(template_sum_squared_error == min(template_sum_squared_error),1);
        
        
        %% Plot :: Audio
        
        figure(audio_fig)
        subplot( n_trial_types,2, 2 * (i_trial_type-1) + 1)
        imagesc(...
            (1:size(peri_stimulus_audio{i_trial_type},2))/sample_rate + psth_flanking_range(1),...
            1:size(peri_stimulus_audio{i_trial_type},1),...
            peri_stimulus_audio{i_trial_type})
        caxis(ylims_audio); colorbar
        %if mod(trial_type_ind,n_trial_type_vs_outcomes) == 1,
        ylabel(trail_type_outcome_yaxis)
        %end
        xlim(psth_flanking_xlim)
        title(i_trial_type)
        
        subplot( n_trial_types, 2, 2 * (i_trial_type-1) + 2)
        imagesc(...
            (1:size(peri_feedback_audio{i_trial_type},2))/sample_rate + psth_flanking_range(1),...
            1:size(peri_feedback_audio{i_trial_type},1),...
            peri_feedback_audio{i_trial_type})
        caxis(ylims_audio); colorbar
        colormap red_white_blue
        xlim(psth_flanking_xlim)
        title(i_trial_type)
        %%
        
        if ~do_all_plots
            continue
        end
        
        
        %% PLOT :: LICKS :: Raster of each trial
        
        %if ~isempty(lick_times)
        figure(lick_raster_fig)
        
        % Smoothing kernel
        lick_rate_ylim = [0 10];
        sigma= 1000/downsample_int; X = -3*sigma:3*sigma;
        kernel = 1/sqrt(2*pi)*sigma * exp(-0.5*X.^2/(sigma^2));
        %kernel = ones(1,sample_rate/5);
        
        timepoint_lick_data = { this_peri_stimulus_licks_all_trials, this_peri_feedback_licks_all_trials};
        timepoint_legend = {'STIMULUS','FEEDBACK'};
        for timepoint_ind = 1:2
            
            subplot( n_trial_types,2, 2 * (i_trial_type-1) + timepoint_ind); cla; hold on
            
            this_timepoint_data = timepoint_lick_data{timepoint_ind};
            [lick_trials,lick_times] = find(this_timepoint_data );
            
            if ~isempty(lick_trials)
                convolved_lick_rate = zeros(size(this_timepoint_data ));
                for i = 1:size(this_timepoint_data ,1)
                    convolved_lick_rate(i,:) = conv(full(this_timepoint_data (i,:)),kernel,'same');
                end
                convolved_lick_rate = mean(convolved_lick_rate);
                [hAx,hLine1,hLine2] = plotyy( ...
                    lick_times/sample_rate  + psth_flanking_range(1), lick_trials,...
                    psth_flanking_times, convolved_lick_rate*(sample_rate/sum(kernel))  );
                set(hAx,{'ycolor'},{'k';'r'})
                % First axis
                set(hLine1,'LineStyle','none','Color','k','Marker','.','MarkerSize',10);
                ylim([0 n_this_trial_type]+[0.5 -0.5]);
                set(hAx(1),'YTick',unique(round(linspace(1,n_this_trial_type-1,3))))
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
            title([ timepoint_legend{timepoint_ind} ' :: ' this_trial_type_legend ' :: ' num2str(i_trial_type)])
        end
        %%
        
        %% PLOT :: Sample trial
        if ~isempty(this_sample_trial)
            
            figure(sample_fig)
            subplot(1,2,1); hold on
            plot( ...
                psth_flanking_times, ...
                this_peri_stimulus_GCaMP_all_trials( this_sample_trial, :),...
                trial_colors{i_trial_type},'LineWidth',2,'markerfacecolor',trial_colors{i_trial_type}(1) );
            title({'Aligned to Stimulus',sprintf('Trial %.0f', this_sample_trial)})
            %ylabel(trial_type_outcome_yaxis)
            ylim(ylims_GCaMP); vline(0,'k')
            xlim(psth_flanking_xlim)
            subplot(1,2,2); hold on
            plot( ...
                psth_flanking_times, ...
                this_peri_feedback_GCaMP_all_trials( this_sample_trial, :),...
                trial_colors{i_trial_type},'LineWidth',2,'markerfacecolor',trial_colors{i_trial_type}(1) );
            title({'Aligned to Stimulus',sprintf('Trial %.0f', this_sample_trial)})
            %ylabel(sprintf('Outcome Type %.0f',trial_group_ind))
            ylim(ylims_GCaMP); vline(0,'k')
            xlim(psth_flanking_xlim);
            
        end


        
        %% PLOT Spectrogram of each trial type
        
        periodogram_fn = @(x) periodogram(nanmean(x,1),...
            hamming(size(x,2)),...
            2.^nextpow2(size(x,2)),sample_rate);
        
        figure(spectrogram_fig);
        subplot( n_trial_types,2, 2 * (i_trial_type-1) + 1); cla; hold on
        % Plot raw in red.
        [pxx,f] = periodogram_fn(this_peri_stimulus_GCaMP_all_trials);
        plot(f,10*log10(pxx),'r')
        % Plot filtered in blue.
        tmp_GCaMP_nansafe = this_peri_stimulus_GCaMP_processed_all_trials;
        tmp_GCaMP_nansafe(isnan(tmp_GCaMP_nansafe )) = nanmean(tmp_GCaMP_nansafe(:));
        periodogram_fn(tmp_GCaMP_nansafe); % Filtered
        h = ylabel(trail_type_outcome_yaxis);
        set(gca,'XScale','log'); legend({'Raw','Filtered'})
        %set(gca,'XTickLabel',arrayfun(@(x) sprintf('%.0f',x*1000),get(gca,'XTick'),'UniformOutput',false))
        
        subplot( n_trial_types, 2, 2 * (i_trial_type-1) + 2); cla; hold on
        [pxx,f] = periodogram_fn(this_peri_feedback_GCaMP_all_trials); % Unfiltered
        plot(f,10*log10(pxx),'r')
        tmp_GCaMP_nansafe = this_peri_feedback_GCaMP_processed_all_trials;
        tmp_GCaMP_nansafe(isnan(tmp_GCaMP_nansafe )) = nanmean(tmp_GCaMP_nansafe(:));
        periodogram_fn(tmp_GCaMP_nansafe); % Filtered
        set(gca,'XScale','log')
        ylabel(''); legend({'Raw','Filtered'})
        %set(gca,'XTickLabel',arrayfun(@(x) sprintf('%.0f',x*),get(gca,'XTick'),'UniformOutput',false))
        
        %% PLOT :: Raster of each trial
        figure(all_trials_fig)
        
        subplot( n_trial_types,2, 2 * (i_trial_type-1) + 1)
        imagesc(...
            psth_flanking_range(1) + (1:size(this_peri_stimulus_GCaMP_processed_all_trials,2))/sample_rate,...
            1:size(this_peri_stimulus_GCaMP_processed_all_trials,1),...
            this_peri_stimulus_GCaMP_processed_all_trials); % Filtered!
        caxis(ylims_GCaMP); colorbar
        ylabel(trail_type_outcome_yaxis)
        xlim(psth_flanking_xlim)
        
        subplot( n_trial_types, 2, 2 * (i_trial_type-1) + 2)
        imagesc( ...
            psth_flanking_range(1) + (1:size(this_peri_feedback_GCaMP_processed_all_trials,2))/sample_rate,...
            1:size(this_peri_feedback_GCaMP_processed_all_trials,1),...
            this_peri_feedback_GCaMP_processed_all_trials)
        caxis(ylims_GCaMP); colorbar
        colormap red_white_blue
        xlim(psth_flanking_xlim);
        
        %% PLOT :: peri-stimulus / peri-feedback AVERAGE
        figure(average_fig)
        aligned_to = {'Stimulus','Feedback'};
        vectors_to_plot = {peri_stimulus_GCaMP_processed_mean,peri_feedback_GCaMP_processed_mean};
        for i_subplot = 1:2
            subplot( n_trial_types,2, 2 * (i_trial_type-1) + i_subplot); hold on
            h = shadedErrorBar( ...
                psth_flanking_times(1:length(vectors_to_plot{i_subplot})), ...
                vectors_to_plot{i_subplot},      ...
                peri_stimulus_GCaMP_processed_err,       ...
                {trial_colors{i_trial_type},'LineWidth',2,'markerfacecolor',trial_colors{i_trial_type}(1)} );
            alpha(h.patch,patch_alpha); delete(h.edge)
            title({'Aligned to Stimulus',sprintf('%.0f trials', n_this_trial_type)})
            
            %ylabel(trial_type_outcome_yaxis)
            %ylim(ylims_GCaMP);
            vline(0,'k');  xlim(psth_flanking_xlim);  ylim(ylims_GCaMP)
        end
        
        
        
        
        %% PLOT :: peri-stimulus average combined figure
        
        figure(combined_fig)
        subplot( 1, 2, 1); hold on;
        h = shadedErrorBar( ...
            psth_flanking_times(1:length(peri_stimulus_GCaMP_processed_mean)), ...
            peri_stimulus_GCaMP_processed_mean,      ...
            peri_stimulus_GCaMP_processed_err,       ...
            {trial_colors{i_trial_type},'LineWidth',2,'markerfacecolor',trial_colors{i_trial_type}(1)} );
        alpha(h.patch,patch_alpha); delete(h.edge)
        title({'Aligned to Stimulus'})
        %ylabel(trial_type_outcome_yaxis)
        ylim(ylims_GCaMP); vline(0,'k')
        
        % Plotting peri-feedback period
        subplot( 1, 2, 2); hold on
        h = shadedErrorBar( ...
            psth_flanking_times(1:length(peri_feedback_GCaMP_processed_mean)),...
            peri_feedback_GCaMP_processed_mean,   ...
            peri_feedback_GCaMP_processed_err,    ...
            {trial_colors{i_trial_type},'LineWidth',2,'markerfacecolor',trial_colors{i_trial_type}(1)} );
        alpha(h.patch,patch_alpha); delete(h.edge)
        
        title({'Aligned to Feedback'})
        %ylabel(sprintf('Outcome Type %.0f',trial_group_ind))
        ylim(ylims_GCaMP); vline(0,'k')
        
        %% Targeted comparisons
        figure(comparison_full_fig)
        
        subplot(3,2,2*(i_trial_outcome-1)+1); hold on
        %         plot(psth_flanking_times(1:length(peri_stimulus_GCaMP_processed_mean)),...
        %             peri_stimulus_GCaMP_processed_mean,  ...
        %             trial_colors{trial_type_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_type_ind})
        h = shadedErrorBar( ...
            psth_flanking_times(1:length(peri_stimulus_GCaMP_processed_mean)), ...
            peri_stimulus_GCaMP_processed_mean,      ...
            peri_stimulus_GCaMP_processed_err,       ...
            {trial_colors{i_trial_type},'LineWidth',2,'markerfacecolor',trial_colors{i_trial_type}(1)} );
        alpha(h.patch,patch_alpha); delete(h.edge)
        ylabel(trial_outcome_legends{i_trial_outcome})
        vline(0); ylim(ylims_GCaMP)
        
        subplot(3,2,2*(i_trial_outcome-1)+2); hold on
        %plot( psth_flanking_times(1:length(peri_feedback_GCaMP_processed_mean)),peri_feedback_GCaMP_processed_mean,  trial_colors{trial_type_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_type_ind})
        h = shadedErrorBar( ...
            psth_flanking_times(1:length(peri_feedback_GCaMP_processed_mean)),...
            peri_feedback_GCaMP_processed_mean,   ...
            peri_feedback_GCaMP_processed_err,    ...
            {trial_colors{i_trial_type},'LineWidth',2,'markerfacecolor',trial_colors{i_trial_type}(1)} );
        alpha(h.patch,patch_alpha); delete(h.edge)
        ylabel(trial_outcome_legends{i_trial_outcome})
        vline(0); ylim(ylims_GCaMP)
        
        
        
        %% Pairwise comparisons
        if do_all_pairwise_plots
            per
            subplot_row_inds = (i_trial_type-1) * n_trial_types + [1:n_trial_types];
            subplot_col_inds = i_trial_type + [0:n_trial_types:(n_trial_types * n_trial_types-1)];
            subplot_all_inds = unique([subplot_row_inds subplot_col_inds]);
            
            figure(pairwise_stimulus_fig)
            for subplot_ind = subplot_all_inds
                subplot( n_trial_types, n_trial_types , subplot_ind); hold on;
                %             shadedErrorBar( ...
                %                 psth_flanking_times, ...
                %                 peri_stimulus_mean,      ...
                %                 peri_stimulus_err,       ...
                %                 {trial_colors{trial_group_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_group_ind}} );
                
                plot(psth_flanking_times(1:length(peri_stimulus_GCaMP_processed_mean)),...
                    peri_stimulus_GCaMP_processed_mean,  ...
                    trial_colors{i_trial_type},...
                    'LineWidth',2,'markerfacecolor',trial_colors{i_trial_type})
                Outcome_Label = sprintf('%s (%.0f)',this_trial_type_legend ,i_trial_type);
                ylim(ylims_GCaMP); vline(0,'k')
                xlim(psth_flanking_xlim)
                
                if subplot_ind <= n_trial_types
                    title(sprintf('%s (%.0f)', Outcome_Label, peri_feedback_n))
                    set(gca,'YTick',[])
                    set(gca,'XTick',[])
                elseif mod(subplot_ind,2) == 1 && (floor(subplot_ind/n_trial_types) == mod(subplot_ind,n_trial_types))
                    axis tight
                else
                    axis off
                end
            end
            
            figure(pairwise_feedback_fig)
            for subplot_ind = subplot_all_inds
                
                % Plotting peri-feedback period
                subplot( n_trial_types, n_trial_types , subplot_ind); hold on;
                %             shadedErrorBar( ...
                %                 psth_flanking_times,...
                %                 peri_feedback_mean,   ...
                %                 peri_feedback_err,    ...
                %                 {trial_colors{trial_group_ind},'LineWidth',2,'markerfacecolor',trial_colors{trial_group_ind}} );
                plot(psth_flanking_times(1:length(peri_feedback_mean)),...
                    peri_feedback_mean,  ...
                    trial_colors{i_trial_type},...
                    'LineWidth',2,'markerfacecolor',trial_colors{i_trial_type})
                
                %ylabel(sprintf('Outcome Type %.0f',trial_group_ind))
                ylim(ylims_GCaMP); vline(0,'k')
                xlim(psth_flanking_xlim)
                
                if subplot_ind <= n_trial_types
                    title(sprintf('%s (%.0f)', Outcome_Label, peri_feedback_n))
                    set(gca,'YTick',[])
                    set(gca,'XTick',[])
                elseif  mod(subplot_ind,2) == 1 && ...
                        (floor(subplot_ind/n_trial_types) == mod(subplot_ind,n_trial_types))
                    axis tight
                    % Leave these axes alone.
                else
                    % Turn off these axes for readability.
                    axis off
                end
                
            end
        end
        
    end % TRIAL TYPE (PLOTTING)
end % Files to loop over.

%tilefigs
disp('Done')
%beep;pause(0.5); beep
