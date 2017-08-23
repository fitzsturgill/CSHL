function [events,epochs] = defineEventsEpochs_Dual2AFCv3
%DEFINEEVENTSEPOCHS_GONOGO   Define events and epochs for spike extraction.
%   [EVENTS,EPOCHS] = DEFINEEVENTSEPOCHS_GONOGO defines events and epochs
%   for spike extraction. 
%
%   EVENTS is a Nx4 cell array with columns corresponding to EventLabel,
%   EventTrigger1, EventTrigger2, Window. EventLabel is the name for
%   referencing the event. EventTrigger1 and EventTrigger2 are names of
%   TrialEvent variables (e.g. 'LeftPortIn'). For fixed windows, the two
%   events are the same; for variable windows, they correspond to the start
%   and end events. Window specifies time offsets relative to the events;
%   e.g. events(1,:) = {'OdorValveOn','OdorValveOn','OdorValveOn',[-3 3]};
%
%   EPOCH is a Nx4 cell array with columns corresponding to  EpochLabel, 
%   ReferenceEvent, Window, RealWindow. EventLabel is the name for 
%   referencing the epoch. ReferenceEvent should match an EventLabel in 
%   EVENTS (used for calculating the epoch rates). RealWindow is currently
%   not implemented (allocated for later versions).
%
%   DEFINEEVENTSEPOCHS_GONOGO defines events and epochs for auditory
%   go-nogo task.
%
%   See also MAKETRIALEVENTS2_GONOGO and DEFINEEVENTSEPOCHS_DEFAULT.

%   Edit log: BH 7/6/12 PM 7/03/14

% Define events and epochs
%              EventLabel       EventTrigger1      EventTrigger2      Window
i = 1;
events(i,:) = {'StimulusOnset',    'StimulusOnset',      'StimulusOnset',      [-6 6]};    i = i + 1;
events(i,:) = {'StimulusOffset',   'StimulusOffset',     'StimulusOffset',     [-6 6]};    i = i + 1;
events(i,:) = {'ResponseStart',    'ResponseStart',      'ResponseStart',      [-6 10]};    i = i + 1;
events(i,:) = {'ResponseEnd',      'ResponseEnd',        'ResponseEnd',          [-10 6]};    i = i + 1;




% Variable events
events(i,:) = {'StimulusSampling',    'StimulusOnset',        'StimulusOffset',     [-6 6]};    i = i + 1;
events(i,:) = {'MovementTime',                'StimulusOffset',       'ResponseStart',     [-6 6]};    i = i + 1;
events(i,:) = {'WaitingRewardTime',                 'ResponseStart',        'ResponseEnd',     [-6 6]};    i = i + 1;


% Define epochs for rate calculations
%               EpochLabel      ReferenceEvent      FixedWindow       RealWindow
i = 1;

% key main epochs
epochs(i,:) = {'StimulusEndRate',    'StimulusOffset',       [-0.25 0.0],     'StimulusSamplingDuration'};    i = i + 1;
epochs(i,:) = {'InitialResponse1',    'ResponseStart',       [0 0.6],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'InitialResponse2',    'ResponseStart',       [0 1],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'InitialResponse3',    'ResponseStart',       [0 4],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'InitialResponse4',    'ResponseStart',       [0 0.35],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'EndWaiting',    'ResponseEnd',       [-0.75 0],       'WaitingRewardTime'};    i = i + 1;
epochs(i,:) = {'RewardRate',    'ResponseEnd',       [0 1],       'WaitingRewardTime'};    i = i + 1;
epochs(i,:) = {'StimulusEndRate',    'StimulusOffset',       [-0.25 0.2],     'StimulusSamplingDuration'};    i = i + 1;

% Variables epochs

epochs(i,:) = {'StimulusResponse','StimulusSampling',[NaN NaN],'NaN'};i = i + 1;
epochs(i,:) = {'MovementResponse','MovementTime',[NaN NaN],'NaN'};i = i + 1;
epochs(i,:) = {'FullWaitingTime','WaitingRewardTime',[NaN NaN],'NaN'};i = i + 1;


end





