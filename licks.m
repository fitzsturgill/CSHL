



file_directory = '/Users/fitz/Documents/KepecsLab/BehaviorData/';

file_name = 'Chat_GCaMP_pink_SO_Training_NIDAQ_Sep22_2015_Session1.mat';
trials_to_skip = [];
trial_outcomes = {};
notch_frequency = NaN;



% file_names{end+1} = 'Chat_GCaMP_pink_SO_Training_NIDAQ_Sep18_2015_Session1.mat';

% notch_frequency = 13; % Some 13hz for Sep02/03, but stronger 20hz for Sep10.




trial_outcomes        = {  [0 4],[1 2], [3 5]};
trial_outcome_fields  = { 'Punish',     'Reward', 'PostUS'      };
trial_outcome_legends = { 'Punishment', 'Reward', 'Omission'    };
trial_colors          = { 'm', 'r' ,  'g','c',      'b','k'           };

% Trial types are 1 (predominantly reward) and 0 (predominantly punishment)
trial_types        = [1 2]; % <-- drawn from data below, but this is usually what
trial_type_legends = {'E(Reward)','E(Punish)'};

lick_raster_fig = namedFigure(['Combined_lick_raster' 'WTF_filter_string']); clf



%% Analysis parameters

psth_flanking_range   = [-2.0 2.1]; % max pre is ~2s because of stimulus presentation.  Can be changed though
psth_flanking_xlim   = [-1.5 2];


%

sample_rate = 10000;

% Y range for all figures
ylims_all = [-2 2];



% Load file

load(fullfile(file_directory, file_name));

fprintf('\nLoaded file %s\n',file_name)

% Trial params
n_trials = length(SessionData.TrialTypes);
n_channels = size(SessionData.NidaqData{1},2);
%     fluorescent_channel = 0 +1; % For Ai0
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




% Loop across trial types
trial_num = 1;
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
%     trial_group_indices  = trial_group_indices(trial_group_indices < trial_too_short_cutoff);
%     fprintf('    %03.0f trials remaining .\n',length(trial_group_indices))
% 
%     fprintf('  trimming with trials_to_skip and trials_to_skip_for_all_fn\n')
%     trial_group_indices = setdiff(trial_group_indices,trials_to_skip);
%     trial_group_indices = setdiff(trial_group_indices,trials_to_skip_for_all_fn(n_trials));
%     trial_group_indices = round(trial_group_indices);
%     fprintf('    %03.0f trials remaining.\n',length(trial_group_indices))

    if length( trial_group_indices ) < 2
        continue
    end



    for trial = trial_group_indices,


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
%         feedback_window = feedback_start_time + psth_flanking_range;

        this_peri_stim_licks     = all_licks( all_licks>stimulus_window(1) & all_licks < stimulus_window(2));
%         this_peri_feedback_licks = all_licks( all_licks>feedback_window (1) & all_licks < feedback_window (2));

        % Convert to samples after beginning of PSTH
        this_peri_stim_licks_psth_samples = ceil(sample_rate*(this_peri_stim_licks-(stimulus_start_time+psth_flanking_range(1))));
%         this_peri_feedback_licks_psth_samples = ceil(sample_rate*(this_peri_feedback_licks-(feedback_start_time+psth_flanking_range(1))));

        % Make a big (sparse!) matrix for this.
        peri_stim_licks{trial_type_outcome_ind}(trial_num, this_peri_stim_licks_psth_samples ) = 1;
%         peri_feedback_licks{trial_type_outcome_ind}(trial_num,this_peri_feedback_licks_psth_samples) = 1;


        %% Advance trial_num
        trial_num = trial_num + 1;

    end % Iterating over trial_group index












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

    timepoint_lick_data = {peri_stim_licks{trial_type_outcome_ind}}; %,peri_feedback_licks{trial_type_outcome_ind}};
    timepoint_legend = {'STIMULUS'} %,'FEEDBACK'}
    for timepoint_ind = 1:1 %:2

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
end
 



