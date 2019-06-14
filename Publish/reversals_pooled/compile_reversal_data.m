%{
script to load and compile data for pooling reversals

you need something like this to prepare prior to calling this script:
DB = dbLoadExperiment('reversals_noPunish_publish');
savepath = fullfile(DB.path, 'pooled', filesep);
smoothWindow = 1;
saveOn = 1;

%%
exp.value = {'DC_44'  'DC_46'  'DC_47' 'DC_53'  'DC_54'  'DC_56'}; % exclude DC_51
exp.valence = {'DC_17'  'DC_20'  'DC_35'  'DC_36'  'DC_37'  'DC_40'};
exp.all = [exp.value exp.valence];
expType = 'all';


%}

%% create structure containing all reversals
AR = struct('csPlus', [], 'csMinus', [], 'csPlusReward', [], 'thirdOdor', []); % AR = all reversals
for si = 1:length(exp.(expType))    
    animal = exp.(expType){si};
%     if strcmp(animal, 'DC_51')
%         continue;
%     end
    load(fullfile(DB.path, 'pooled', ['RE_' animal '.mat']));
    for group = fieldnames(AR)'
        sgroup = group{:};
        for field = fieldnames(RE.(sgroup))'
            sfield = field{:};          
            if any(strcmp(sfield, {'trialsBefore', 'trialsAfter'})) % these are special fields that don't contain before and after data
                continue
            end
            if si == 1
                AR.(sgroup).(sfield).before = RE.(sgroup).(sfield).before;
                AR.(sgroup).(sfield).after = RE.(sgroup).(sfield).after;
            else
                AR.(sgroup).(sfield).before = expandVertCat(AR.(sgroup).(sfield).before, RE.(sgroup).(sfield).before, 'right');
                AR.(sgroup).(sfield).after = expandVertCat(AR.(sgroup).(sfield).after, RE.(sgroup).(sfield).after, 'left');
            end
        end
    end
end

firstReversals = cellfun(@(x,y) ~strcmp(x(1:5), y(1:5)), AR.csPlus.filename.before(1:end-1,end), AR.csPlus.filename.before(2:end,end));
firstReversals = [false; firstReversals];
nReversals = size(AR.csPlus.globalTrialNumber.after, 1);
% reversal #s
revNumber = ones(nReversals,1);
mouseNumber = ones(nReversals,1);

thisRev = 1;
thisMouse = 1;
for counter = 2:nReversals
    if strcmp(AR.csPlus.filename.before{counter - 1,end}(1:5), AR.csPlus.filename.before{counter,end}(1:5))
        thisRev = thisRev + 1;
    else
        thisRev = 1;
        thisMouse = thisMouse + 1;
    end
    revNumber(counter) = thisRev;
    mouseNumber(counter) = thisMouse;
end


%% compile data
fieldsToCompile = {...
    'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'phPeakPercentile_cs_ch1', 'phPeakPercentile_cs_ch2', 'csLicksROC', 'licks_cs', 'pupil_cs', 'pupil_csBaselined' 'whisk_cs', 'wheel_baseline',...
    'phPeakMean_us_ch1_deconv', 'phPeakMean_us_ch2_deconv', 'phPeakPercentile_us_ch1_deconv', 'phPeakPercentile_us_ch2_deconv', 'phBaseline_ch1', 'phBaseline_ch2', 'phPeakMean_us_ch1', 'phPeakMean_us_ch2'};

newCsPlus = struct();
newCsMinus = struct();
alwaysCsPlus = struct();
alwaysCsPlusReward = struct();
odor3 = struct();

for counter = 1:length(fieldsToCompile)
    field = fieldsToCompile{counter};
    newCsPlus.(field) = smoothdata([AR.csMinus.(field).before AR.csPlus.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
    newCsMinus.(field) = smoothdata([AR.csPlus.(field).before AR.csMinus.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
    alwaysCsPlus.(field) = smoothdata([AR.csPlus.(field).before AR.csPlus.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
    alwaysCsPlusReward.(field) = smoothdata([AR.csPlusReward.(field).before AR.csPlusReward.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
    odor3.(field) = smoothdata([AR.thirdOdor.(field).before AR.thirdOdor.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
end

% newCsPlus.trialNumber =  
newCsPlus.trialNumber = (1:size(newCsPlus.licks_cs, 2)) - size(AR.csMinus.licks_cs.before, 2);
newCsPlus.firstRevTrial = size(AR.csMinus.licks_cs.before, 2) + 1; 
newCsMinus.trialNumber = (1:size(newCsMinus.licks_cs, 2)) - size(AR.csPlus.licks_cs.before, 2);  
newCsMinus.firstRevTrial = size(AR.csPlus.licks_cs.before, 2) + 1; 
alwaysCsPlus.trialNumber = (1:size(alwaysCsPlus.licks_cs, 2)) - size(AR.csPlus.licks_cs.before, 2); 
alwaysCsPlus.firstRevTrial = size(AR.csPlus.licks_cs.before, 2) + 1; 
alwaysCsPlusReward.trialNumber = (1:size(alwaysCsPlusReward.licks_cs, 2)) - size(AR.csPlusReward.licks_cs.before, 2); 
alwaysCsPlusReward.firstRevTrial = size(AR.csPlusReward.licks_cs.before, 2) + 1; 
odor3_trialNumber = (1:size(odor3.licks_cs, 2)) - size(AR.thirdOdor.licks_cs.before, 2); 

oldCsPlus_trialNumber = -(size(AR.csPlus.licks_cs.before, 2) - 1) : 0;
oldCsMinus_trialNumber = -(size(AR.csMinus.licks_cs.before, 2) - 1) : 0;

%%
% quality control- calculate auROC and dPrime for relevent comparisons
trialWindow = [-20 60];%?

% initialize
comp = {'licks_cs', 'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'pupil_csBaselined', 'whisk_csBaselined'}; % comparisons
all_ways = struct(...
    'before', zeros(nReversals, 1),... % compare cs- and cs+ before reversal
    'after', zeros(nReversals, 1),... % compare cs- and cs+ after reversal
    'acq', zeros(nReversals, 1),... % acquisition: compare cs- before and cs+ after reversal
    'ext', zeros(nReversals, 1)...  % extinction: compare cs+ before and cs- after reversal
    );
% clear dPrime auROC
for field = comp
    dPrime.(field{:}) = all_ways;
    auROC.(field{:}) = all_ways;
end

% calculate metrics
for field = comp
    for rev = 1:nReversals
        % before
        dPrime.(field{:}).before(rev) = dPrime_SNR(AR.csPlus.(field{:}).before(rev,trialWindow(1) + end + 1:end), AR.csMinus.(field{:}).before(rev,trialWindow(1) + end + 1:end));
        auROC.(field{:}).before(rev) = rocarea(stripNaNs(AR.csPlus.(field{:}).before(rev,trialWindow(1) + end + 1:end)), stripNaNs(AR.csMinus.(field{:}).before(rev,trialWindow(1) + end + 1:end)), 'scale');
        % after
        dPrime.(field{:}).after(rev) = dPrime_SNR(AR.csPlus.(field{:}).after(rev,1:trialWindow(2)), AR.csMinus.(field{:}).after(rev,1:trialWindow(2)));
        auROC.(field{:}).after(rev) = rocarea(stripNaNs(AR.csPlus.(field{:}).after(rev,1:trialWindow(2))), stripNaNs(AR.csMinus.(field{:}).after(rev,1:trialWindow(2))), 'scale');        
        % acquisition
        dPrime.(field{:}).acq(rev) = dPrime_SNR(AR.csPlus.(field{:}).after(rev,1:trialWindow(2)), AR.csMinus.(field{:}).before(rev,trialWindow(1) + end + 1:end));
        auROC.(field{:}).acq(rev) = rocarea(stripNaNs(AR.csPlus.(field{:}).after(rev,1:trialWindow(2))), stripNaNs(AR.csMinus.(field{:}).before(rev,trialWindow(1) + end + 1:end)), 'scale');                
        % extinction
        dPrime.(field{:}).ext(rev) = dPrime_SNR(AR.csPlus.(field{:}).before(rev,trialWindow(1) + end + 1:end), AR.csMinus.(field{:}).after(rev,1:trialWindow(2)));
        auROC.(field{:}).ext(rev) = rocarea(stripNaNs(AR.csPlus.(field{:}).before(rev,trialWindow(1) + end + 1:end)), stripNaNs(AR.csMinus.(field{:}).after(rev,1:trialWindow(2))), 'scale');
    end    
end

% quality control % 2
% detect when auROC values first exceed threshold
rocThresh = 0.5;
trialsToCriterion = NaN(nReversals, 1);
% looping is just easier
for counter = 1:nReversals    
    thisRev = AR.csPlus.csLicksROC.after(counter, :);
    thisRev = thisRev > rocThresh;
    nt = find(thisRev, 1);
    if ~isempty(nt)
        trialsToCriterion(counter) = nt;
    end
end

