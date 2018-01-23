%% desktop 
files = {...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_17\', 'RE_DC_17.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_18\', 'RE_DC_18.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_20\', 'RE_DC_20.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_25\', 'RE_DC_25.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_27\', 'RE_DC_27.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_36\', 'RE_DC_36.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_40\', 'RE_DC_40.mat';...
    };


%%
for counter = 1:size(files, 1)
    fileName = files{counter, 2}(1:end-4);
    temp = load(fullfile(files{counter, 1}, files{counter, 2}));
    if counter == 1
        RE = temp.RE;
    else
        RE(counter) = temp.RE;
    end
end

%% 
%% desktop
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SFN_2017\Reversals';
savepath = 'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal';
saveOn = 1;

%% create structure containing all reversals
AR = struct('csPlus', [], 'csMinus', [], 'csPlusReward', [], 'allTrials', []); % AR = all reversals
for group = fieldnames(AR)'
    sgroup = group{:};
    for field = fieldnames(RE(1).(sgroup))'
        sfield = field{:};
        for si = 1:length(RE)
            if sum(strcmp(sfield, {'trialsBefore', 'trialsAfter'})) % these are special fields that don't contain before and after data
                continue
            end
            if si == 1
                AR.(sgroup).(sfield).before = RE(si).(sgroup).(sfield).before;
                AR.(sgroup).(sfield).after = RE(si).(sgroup).(sfield).after;
            else
                AR.(sgroup).(sfield).before = expandVertCat(AR.(sgroup).(sfield).before, RE(si).(sgroup).(sfield).before, 'right');
                AR.(sgroup).(sfield).after = expandVertCat(AR.(sgroup).(sfield).after, RE(si).(sgroup).(sfield).after, 'left');
            end
        end
    end
end
nReversals = size(AR.csPlus.globalTrialNumber.after, 1);