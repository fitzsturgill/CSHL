sample_rate = 10000;
trial_outcome_fields  = { 'Punish',     'Reward', 'PostUS'      };
    % PSTH samples.
psth_flanking_range   = [-2.0 2.1]; % max pre is ~2s because of stimulus presentation.  Can be changed though
psth_flanking_xlim   = [-1.5 2];

psth_flanking_times   = linspace( psth_flanking_range(1),psth_flanking_range(2), sample_rate*diff(psth_flanking_range)+1  );
psth_flanking_samples = psth_flanking_times * sample_rate;
feedback_start_time = SessionData.RawEvents.Trial{trial}.States.(trial_outcome_fields{trial_outcome_ind})(1);


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

