% for paul
DB = dbLoadExperiment('reversals_noPunish_publish');
animal = 'DC_47';
dbLoadAnimal(DB, animal);

savepath = fullfile(DB.path, ['pooled' filesep 'ForPaul']);
ensureDirectory(savepath);

fieldsToCopy = {'PhotometryExpFit', 'Whisk', 'Wheel', 'pupil', 'filename', 'trialNumber', 'trialOutcome', 'trialType', 'sessionIndex',...
    'ITI', 'OdorValveIndex', 'BlockNumber', 'ReinforcementOutcome', 'reject', 'TrialStartTimestamp'};

allData = struct();
for counter = 1:length(fieldsToCopy)
    field = fieldsToCopy{counter};
    allData.(field) = TE.(field);
end
[lickTimes, lickTrials] = extractEventTimesFromTE(TE, 1:length(TE.trialNumber), 'Port1In', 'zeroField', 'PreCsRecording', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
allData.licks.lickTimes = lickTimes;
allData.licks.lickTrials = lickTrials;
save(fullfile(savepath, sprintf('fitzData_%s.mat', animal)), 'allData');







%% FrankenLNL_RewardPunish



animal = 'ACh_13';
filepath = 'Z:\SummaryAnalyses\Franken_LNL_aversive\ACh_13_FrankenLNL_4odors_Apr25_2019_Session1\';
% load(fullfile(filepath, 'TE.mat'));
frankenLNL_conditions;

fieldsToCopy = {'PhotometryExpFit', 'Whisk', 'Wheel', 'pupil', 'filename', 'trialNumber', 'trialOutcome', 'trialType', 'sessionIndex',...
    'ITI', 'Odor2ValveIndex', 'BlockNumber', 'ReinforcementOutcome', 'reject', 'TrialStartTimestamp'};

allData = struct();
for counter = 1:length(fieldsToCopy)
    field = fieldsToCopy{counter};
    allData.(field) = TE.(field);
end
[lickTimes, lickTrials] = extractEventTimesFromTE(TE, 1:length(TE.trialNumber), 'Port1In', 'zeroField', 'PreCsRecording', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
allData.licks.lickTimes = lickTimes;
allData.licks.lickTrials = lickTrials;
save(fullfile(savepath, sprintf('fitzData_%s.mat', animal)), 'allData');



%% 

DB = dbLoadExperiment('FrankenLNL_RewardPunish');

animals = {'ACh_3', 'ACh_7', 'ACh_13'};

animal = animals{3};
dbLoadAnimal(DB, animal);

