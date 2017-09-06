% LNL_odor_v2_pav_rev  reversal analysis script 9/
%%
files = {...
    'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_17\', 'RE_DC_17.mat';...
    'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_20\', 'RE_DC_20.mat';...
    };

for counter = 1:size(files, 1);
    fileName = files{counter, 2}(1:end-4);
    temp = load(fullfile(files{counter, 1}, files{counter, 2}));
    if counter == 1;
        RE = temp.RE;
    else
        RE(counter) = temp.RE;
    end
end

%% 
% find total number of reverals, max post reversal trials, etc.

