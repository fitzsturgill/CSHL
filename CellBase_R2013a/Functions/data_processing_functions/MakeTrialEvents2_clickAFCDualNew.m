function MakeTrialEvents2_clickAFCDualNew(sessionpath)
%MAKETRIALEVENTS2_GONOGO   Synchronize trial events to recording times. 
%	MAKETRIALEVENTS2_GONOGO(SESSIONPATH) loads Neuralynx events and adjusts
%	trial event times (trial-based behavioral data, see
%	SOLO2TRIALEVENTS4_AUDITORY_GONOGO) to the recorded time stamps. This
%	way the neural recordings and behavioral time stamps are in register.
%	Stimulus time TTL pulses are used for synchronization. The synchronized
%	trial events structure is saved under the name 'TrialEvents.mat'. This
%	file becomes the primary store of behavioral data for a particular
%	session; it is retrieved by LOADCB via CELLID2FNAMES. This default
%	file name is one of the preference settings of CellBase - type 
%   getpref('cellbase','session_filename');
%
%   MAKETRIALEVENTS2_GONOGO(SESSIONPATH,'StimNttl',TTL) specifies the TTL
%   channel which serves as the basis for synchronization.
%
%   See also SOLO2TRIALEVENTS4_AUDITORY_GONOGO and LOADCB.

% Parse input arguments

% Load converted Neuralynx event file 
if ~isdir(sessionpath),
    error('Session path is WRONG');
end
try
    load(fullfile(sessionpath, 'Events.mat'));   % load Neuralynx events
catch    
    error('EVENTS file not found; make sure Nlx files have been converted.');
end
try
    load(fullfile(sessionpath, 'EventTTL_task.mat'));   % load Neuralynx events
catch    
    error('EVENTS TTL file not found; make sure Nlx files have been converted.');
end

EventTimestamps=EventTimestamps; % EventTimestamps is the name of the variable in Events.mat

% Load trial events structure
SE_filename = [sessionpath filesep 'TE.mat'];
TE = load(SE_filename);

EventTTL_task=EventTTL;

% Create cell aray to group TTLs adn their timestamps by trial. 
[NlxEvents] = GetTrialStartTimeStampPM2(EventTTL_task, EventTimestamps);

% Set TTl alignement state at 2 Waitingfor initial poke) Index 4 is DeliverStimulus1

idx=2; 

AlignedNlxEvents=cell(3,length(NlxEvents));

for i=1:length(NlxEvents)
    AlignedNlxEvents{1,i}=NlxEvents{1,i};
    size(NlxEvents{2,i}(AlignedNlxEvents{1,i}==idx));
    AlignedNlxEvents{2,i}=NlxEvents{2,i}-NlxEvents{2,i}(AlignedNlxEvents{1,i}==idx);
    AlignedNlxEvents{3,i}=NlxEvents{2,i}(AlignedNlxEvents{1,i}==idx);
end


% Synchronization
son = find(EventTTL_task==idx); 

TE2 = TE.TE;
son2 = EventTimestamps(son)/(10^6);   % Trial start time recorded by the recording system (Neuralynx)
ts = TE2.TrialStartTimeStamp;   % Trial start in absolut time recorded by the behavior control system



% Match timestamps - in case of mismatch, try to fix
if ~ismatch(ts,son2)
    % note: obsolete due the introduction of TTL parsing
    son2 = clearttls(son2); % eliminate recorded TTL's within 0.5s from each other - broken TTL pulse
    if ~ismatch(ts,son2)
        son2 = trytomatch(ts,son2);  % try to match time series by shifting
        if ~ismatch(ts,son2)
            son2 = tryinterp(ts,son2); % interpolate missing TTL's or delete superfluous TTL's up to 10 erroneous TTl's
            if ~ismatch(ts,son2)  % TTL matching failure
                error('MakeTrialEvents:TTLmatch','Matching TTLs failed.')
            else
                warning('MakeTrialEvents:TTLmatch','Missing TTL interpolated.')
            end
        else
            warning('MakeTrialEvents:TTLmatch','Shifted TTL series.')
        end
    else
        warning('MakeTrialEvents:TTLmatch','Broken TTLs cleared.')
    end
end

% Eliminate last TTL's recorded in only one system
sto = TE2.TrialStartTimeStamp;
if length(son2) > length(ts)   % time not saved in behavior file (likely reason: autosave was used)
    son2 = son2(1:length(ts));
elseif length(son2) < length(ts)  % time not recorded on Neuralynx (likely reason: recording stopped)
    shinx = 1:length(son2);
    ts = ts(shinx);
    sto = sto(shinx);
    TE2
    shinx
    TE2 = shortenTE(TE2,shinx);
end

TE2.TrialStartAligned = son2;


% Save synchronized 'TrialEvents' file
save([sessionpath filesep 'TrialEvents.mat'],'-struct','TE2')

if ~isempty(TE2.TrialStartTimeStamp),
    save([sessionpath filesep 'TrialEvents.mat'],'-struct','TE2')
else
    error('MakeTrialEvents:noOutput','Synchronization process failed.');
end

if ~isempty(AlignedNlxEvents),
    save([sessionpath filesep 'AlignedNlxEvents.mat'],'AlignedNlxEvents')

end

% -------------------------------------------------------------------------
function I = ismatch(ts,son2)

% Check if the two time series match notwithstanding a constant drift
clen = min(length(ts),length(son2));
I = abs(max(diff(ts(1:clen)-son2(1:clen)))) < 3;  % the difference between the timestamps on 2 systems may have a constant drift, but it's derivative should still be ~0

% note: abs o max is OK, the derivative is usually a small neg. number due
% to drift of the timestamps; max o abs would require a higher tolerance
% taking the drift into account (if 2 event time stamps are far, the drift
% between them can be large)

% -------------------------------------------------------------------------
function son2 = tryinterp(ts,son2)

% Interpolate missing TTL's or delete superfluous TTL's up to 10 erroneous TTl's
for k = 1:10
    if ~ismatch(ts,son2)
        son3 = son2 - son2(1) + ts(1);
        adt = diff(ts(1:min(length(ts),length(son2)))-son2(1:min(length(ts),length(son2))));
        badinx = find(abs(adt)>0.1,1,'first') + 1;  % find problematic index
        if adt(badinx-1) < 0    % interploate
            ins = ts(badinx) - linterp([ts(badinx-1) ts(badinx+1)],[ts(badinx-1)-son3(badinx-1) ts(badinx+1)-son3(badinx)],ts(badinx));
            son2 = [son2(1:badinx-1) ins+son2(1)-ts(1) son2(badinx:end)];
        else
%             ins = son3(badinx) - linterp([son3(badinx-1) son3(badinx+1)],[son3(badinx-1)-ts(badinx-1) son3(badinx+1)-ts(badinx)],son3(badinx));
%             ts = [ts(1:badinx-1) ins ts(badinx:end)];
            son2(badinx) = [];   % delete
        end
    end
end

% -------------------------------------------------------------------------
function son2 = trytomatch(ts,son2)

% Try to match time series by shifting
len = length(son2) - 15;
minx = nan(1,len);
for k = 1:len
    minx(k) = max(diff(ts(1:15)-son2(k:k+14)));  % calculate difference in the function of shift
end
mn = min(abs(minx));
minx2 = find(abs(minx)==mn);
minx2 = minx2(1);   % find minimal difference = optimal shift
son2 = son2(minx2:min(minx2+length(ts)-1,length(son2)));

% -------------------------------------------------------------------------
function son2 = clearttls(son2)

% Eliminate recorded TTL's within 0.5s from each other
inx = [];
for k = 1:length(son2)-1
    s1 = son2(k);
    s2 = son2(k+1);
    if s2 - s1 < 0.5
        inx = [inx k+1]; %#ok<AGROW>
    end
end
son2(inx) = [];

% -------------------------------------------------------------------------
function TE2 = shortenTE(TE2,shinx)

% Eliminate behavioral trials
fnm = fieldnames(TE2);
for k = 1:length(fieldnames(TE2))
    TE2.(fnm{k}) = TE2.(fnm{k})(shinx);
end



