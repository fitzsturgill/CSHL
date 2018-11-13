% for paul
DB = dbLoadExperiment('reversals_noPunish_publish');
animal = 'DC_44';
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




