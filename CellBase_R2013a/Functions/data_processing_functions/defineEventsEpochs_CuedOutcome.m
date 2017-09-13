function [events,epochs] = defineEventsEpochs_CuedOutcome
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


% Define events and epochs
%              EventLabel       EventTrigger1      EventTrigger2      Window
i = 1;
events(i,:) = {'Cue_start',    'Cue_start',      'Cue_start',      [-4 7]};    i = i + 1; % aligned to cue onset
events(i,:) = {'Us_start',    'Us_start',      'Us_start',      [-7 4]};    i = i + 1;


% Variable events


% Define epochs for rate calculations
%               EpochLabel      ReferenceEvent      FixedWindow       RealWindow
i = 1;
epochs(i,:) = {'Baseline',    'Cue_start',       [-4 0],        ''};    i = i + 1;
epochs(i,:) = {'Phasic',    'Cue_start',       [0.25 1.25],        ''};    i = i + 1;
epochs(i,:) = {'Sustained',    'Cue_start',       [1.25 3],        ''};    i = i + 1;
epochs(i,:) = {'Outcome',    'Us_start',       [0 1.5],        ''};    i = i + 1;
