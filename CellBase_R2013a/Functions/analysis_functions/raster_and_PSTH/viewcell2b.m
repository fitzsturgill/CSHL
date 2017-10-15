function viewcell2b(cellid,varargin)
%VIEWCELL2B   Raster plot and PSTH.
%   VIEWCELL2B(CELLID,'TRIGGERNAME',TRIGEVENT,'SORTEVENT',SEVENT,'EVENTTYPE
%   ',EVTYPE) plots a light triggered or event triggered raster and PSTH
%   for a given cell. TRIGEVENT is used as trigger and trials are sorted
%   according to SEVENT. 'EVENTTYPE' can be 'STIM' or 'BEHAV'. Optional
%   input argument pairs:
%       'ShowEvets', events indicated on the raster plots
%       'ShowEventsColors', specifies colors for the events shown
%       'FigureNum', specifies the figure handle to plot on
%       'window', window size for the plots
%       'dt', time resoultion
%       'sigma', determines the smoothing kernel for the smoothed PSTH (see
%           SMOOTHED_PSTH)
%       'isadaptive', default = false - if true, adaptive PSTH algorithm is
%           used (see APSTH and BINRASTER2APSTH)
%       'PSTHstd', 'on' or 'off', controls whether SD is plotted on the PSTH
%           panel
%       'Partitions', e.g. 'all', '#ResponseType', can be used to partition
%           the trials according to different variables, multiple rasters 
%           are plotted
%       'EventMarkerWidth', specifies marker size for events shown
%           'PlotZeroLine', 'on' or 'off', controls whether a line
%           indicating zero appears on the plots
%
%   For example calls, see QUICKANALISYS2.
%
%   See also PREALIGNSPIKES, STIMES2BINRASTER, BINRASTER2PSTH,
%   PLOT_RASTER2A.

%   Edit log: BH 6/23/11, 2/10/12

% Input argument check
if nargin < 1
    help viewcell2b
    return
end

% Default arguments
default_args={...
    'window',               [-0.5 1];...
    'dt',                   0.01;...
    'sigma',                0.02;...
    'isadaptive'            false;...
    'FigureNum',            1;...
    'TriggerName',          'PulseOn';...
    'SortEvent',            'PulseOn';...
    'eventtype',             'stim';... % 'behav'
    'ShowEvents',           {'PulseOn'};...
    'ShowEventsColors',     {'r'};...
    'Num2Plot',             'all';...
    'PlotDashedEvent',      '';...
    'PlotDashedCondition',  'min';...
    'PSTHPlot',             1;...
    'PSTHlinewidth',        1.5;...
    'DashedLineStyle',      ':';...
    'LastEvents',           '';...
    'Partitions',           'all';...
    'PrintCellID',          'on';...
    'PrintCellIDPos',       'bottom-right';...
    'BurstPSTH'             'off';......  
    };
[g,error] = parse_args(default_args,varargin{:});

% Check if cellid is valid
if validcellid(cellid,{'list'}) ~= 1
    fprintf('%s is not valid.',cellid);
    return
end

%--------------------------------------------------------------------------
% Preprocessing
%--------------------------------------------------------------------------

% Time
margin = g.sigma * 3;     % add an extra margin to the windows
time = g.window(1)-margin:g.dt:g.window(2)+margin;   % time base array

% Load events
switch g.eventtype
    case 'stim'
        TE = loadcb(cellid,'StimEvents');
        SP = loadcb(cellid,'STIMSPIKES');
    case {'event','behav'}
        TE = loadcb(cellid,'TrialEvents'); % 'struct' option removed in 8.6.0.267246 (R2015b)
        SP = loadcb(cellid,'EVENTSPIKES');
end
trigger_pos = findcellstr(SP.events(:,1),g.TriggerName);

% TriggerName mismatch
if trigger_pos == 0
    error('Trigger name not found');
else
    g.TriggerEvent = SP.events{trigger_pos,2};
end
if ~isfield(TE,g.TriggerEvent)
    error('TriggerEvent mismatch: supply correct Events structure')
end

% Spike times
alltrials = 1:size(SP.event_stimes{1},2);
stimes  = SP.event_stimes{trigger_pos}(alltrials);
if strcmp(g.BurstPSTH,'on'),
    stimes = detect_bursts(stimes);
end

% Event windows
if ~iscellstr(g.LastEvents) && (strcmpi(g.LastEvents,'none') || isempty(g.LastEvents))
    window_margin = SP.events{trigger_pos,4};
    ev_windows = SP.event_windows{trigger_pos};
else
    window_margin = [g.window(1)-2*g.dt 0];
    ev_windows = get_last_evtime(TE,g.TriggerEvent,g.LastEvents);
end

%--------------------------------------------------------------------------
% Make the main raster
%--------------------------------------------------------------------------

% Calculate binraster
NUMtrials = length(TE.(g.SortEvent));
if iscell(stimes{1})   % deal with lick-aligned raster
    stimes2 = [stimes{1:end}];
    binraster0 = stimes2binraster(stimes2,time,g.dt);
    binraster = nan(NUMtrials,size(binraster0,2));
%     binraster2 = nan(NUMtrials,size(binraster0,2));
    for k = 1:NUMtrials   % calculate sum of rows for each trial, which will be used for the PSTH
        sind = sum(cellfun(@length,stimes(1:k-1))) + 1;
        eind = sind + length(stimes{k}) - 1;
%         disp([sind eind])
        binraster(k,:) = mean(binraster0(sind:eind,:),1);
%         binraster2(k,:) = sum(stimes2binraster(stimes{k},time,g.dt),1);
    end
else
    binraster = stimes2binraster(stimes,time,g.dt);
end

% Partition trials
[COMPTRIALS, TAGS] = partition_trials(TE,g.Partitions);
vinx = cellfun(@(s)(~isempty(s)),COMPTRIALS);
COMPTRIALS = COMPTRIALS(vinx);
TAGS = TAGS(vinx);
trigev = TE.(g.TriggerEvent);
if ~iscell(trigev)
    valid_trials = find(~isnan(trigev));
else
    valid_trials = find(cellfun(@(s)~isempty(s),trigev));
end

% Calculate PSTH
if ~g.isadaptive
    [psth, spsth, spsth_se] = binraster2psth(binraster,g.dt,g.sigma,COMPTRIALS,valid_trials);
else
    [psth, spsth, spsth_se] = binraster2apsth(binraster,g.dt,g.sigma,COMPTRIALS,valid_trials);   % adaptive PSTH
end
EventTimes = trialevents2relativetime(TE,g.TriggerEvent,g.ShowEvents);

% Sort trials
NUMevents = length(g.SortEvent);
if iscellstr(g.SortEvent)
    sort_var = nan(NUMevents,NUMtrials);
    for iS = 1:NUMevents
        sort_var(iS,:) = TE.(g.SortEvent{iS}) - TE.(g.TriggerEvent);
    end
    sort_var = min(sort_var);
elseif ~isempty(g.SortEvent)
     if ~iscell(TE.(g.TriggerEvent))
        sort_var = TE.(g.SortEvent) - TE.(g.TriggerEvent);
    else
        gte = nan(1,NUMtrials);
        inx = ~cellfun(@isempty,TE.(g.TriggerEvent));
        gte(inx) = cell2mat(cellfun(@(s)s(1),TE.(g.TriggerEvent)(inx),'UniformOutput',false));
        sort_var = TE.(g.SortEvent) - gte;
    end
else
    sort_var = NaN;
end
[mylabels, mycolors, mycolors2,mylinestyle] = makeColorsLabels(@defineLabelsColors_Balazs,TAGS);
XLabel = ['Time - ' g.TriggerEvent];
YLabel = 'Rate (Hz)';

%-----------------------------
% Raster + PSTH
%-----------------------------

% Plot the raster
fhandle0 = plot_raster2a(stimes,time,valid_trials,COMPTRIALS,mylabels,EventTimes,window_margin,ev_windows,sort_var,...
    g,'Colors',mycolors,'Colors2',mycolors2,'NumTrials2Plot',g.Num2Plot);
% fhandle0 = plot_raster2a(stimes,time,valid_trials,COMPTRIALS,mylabels,EventTimes,window_margin,ev_windows,sort_var,...
%     varargin,'Colors',{mycolors},'Colors2',{mycolors2},'NumTrials2Plot',g.Num2Plot);

if isfield(g,'Legend')
    mylabels = g.Legend;
end

% Plot the PSTH
if g.PSTHPlot == 1
    if ~isempty(g.PlotDashedEvent)
        if nansum(TE.TrialStart) ~= 0   % i.e. TrialStart is an actual timestamp
            g.PlotDashedTime = eval(['TE.' g.PlotDashedEvent]);
        else   % if TrialStart is a dummy variable set to 0
            g.PlotDashedTime = SP.event_windows{trigger_pos}(2,valid_trials);
        end
        switch g.PlotDashedCondition
            case 'min'
                g.PlotDashedTime = nanmin(g.PlotDashedTime);
            case 'max'
                g.PlotDashedTime = nanmax(SP.event_windows{trigger_pos}(2,valid_trials));
            case 'mean'
                g.PlotDashedTime = nanmean(SP.event_windows{trigger_pos}(2,valid_trials));
            case 'median'
                g.PlotDashedTime = nanmedian(SP.event_windows{trigger_pos}(2,valid_trials));
        end
    end
    plot_timecourse(time,spsth,spsth_se,g,'FigureNum',fhandle0(end),'Colors',mycolors,'LineStyle',mylinestyle,'Legend',mylabels,'XLabel',XLabel,'YLabel',YLabel);
    axis tight
end
if strcmpi(g.PrintCellID,'on')
    fstamp(cellid,80,g.PrintCellIDPos);
end