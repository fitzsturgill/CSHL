% multilinear regression with licking

% load('Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_37\TE.mat');
load('Z:\SummaryAnalyses\LickNoLick_odor_v2\DC_35\TE.mat');
LNL_conditions;
%%
TE.licksTS = bpEventToTimeSeries(TE, 'Port1In', 'duration', 11);

%%
% Initial Conditions
Fs = 20; % sample rate
channel = 1;
Kw = [-2 2]; % Kernel/weight offset window


% generate weight time offset indices
Kw2 = Kw * Fs;
Kwix = Kw2(1):Kw2(2);
margins = [sum(Kwix < 0) sum(Kwix) > 0]; % to avoid edge effects due to using circ-shifted predictor variables

% Predictor
T = TE.Photometry.data(channel).ZS';
T = T(margins(1) + 1 : end - margins(2), :);
truncSize = size(T);
% T = ([1 1 1 1; 2 2 2 2; 3 3 3 3]);
T = T(:);
T = zscore(T);

% Design matrix
dm = TE.licksTS.data' > 0; % get rid of spurious fast bursts of licks
% dm = ([0:100; 50:150; 100:200])';
DM = zeros(truncSize(1), truncSize(2), length(Kwix));

% shift matrices
for counter = 1:length(Kwix)
    shifted = circshift(dm, Kwix(counter), 1);
    DM(:,:,counter) = shifted(margins(1) + 1 : end - margins(2), :);
end

% reshape
DM = reshape(DM, prod(truncSize), length(Kwix));

% regress
[w wconf R] = regress(T, DM);

% xData for plotting weight kernel
xData = linspace(Kw(1), Kw(2), length(w));
%
ensureFigure('test', 1); boundedline(xData', w, wconf - w);
ensureFigure('test', 1); plot(xData', w, 'b'); hold on; plot(xData', wconf,'c')
R = reshape(R, truncSize(1), truncSize(2));
R = R';
%
ensureFigure('test2', 1); imagesc(R, [-2 2]);
