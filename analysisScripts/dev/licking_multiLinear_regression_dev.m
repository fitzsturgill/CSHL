% multilinear regression with licking

load('Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_37\TE.mat');
LNL_conditions;


%%
% Initial Conditions
Fs = 20; % sample rate
channel = 1;

% Predictor
T = TE.Photometry.data(channel).ZS';

% T = ([1 1 1 1; 2 2 2 2; 3 3 3 3]); 
T = T(:);
T = zscore(T);

% Kernel/weight offset window
Kw = [-2 2];
Kw = Kw * Fs;
% weight time offset indices
Kwix = Kw(1):Kw(2);

% Design matrix
dm = TE.licksTS.data' > 0;
% dm = ([0:100; 50:150; 100:200])'; 
DM = zeros(size(dm,1), size(dm,2), length(Kwix));

% shift matrices
for counter = 1:length(Kwix)
    DM(:,:,counter) = circshift(dm, Kwix(counter), 1);
end

% reshape
DM = reshape(DM, numel(dm), length(Kwix));

% regress
W = regress(T, DM);

%%
ensureFigure('test', 1); plot(W);

disp('good progress but see commit, issue of predictors circ shift edge effects');
