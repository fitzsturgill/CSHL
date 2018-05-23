%% experiment 1-  mimic of my reversal learning task
savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\Modeling\KalmanRW_Reversal_Task';

% initialize using Gershman code, prior variance for weights = 1,
% variance of process noise = 0.01, variance of measurement noise = 1
param = KTD_defparam;
rewardSize = 1;
punishSize = -1;
revFreq = 100; % block size (normal or reversed contingencies alternate)
n = 500; % number of trials
trials = (1:n)';
cueType = randi(2,n,1); % each trial contains either stimulus 1 or 2
X = zeros(n,2); % stimuli matrix nTrials x 2,  
X(cueType == 1, 1) = 1;
X(cueType == 2, 2) = 1;


reversed = floor(trials / revFreq);
reversed = rem(reversed, 2) == 1;

% randomly deliver reward on 80% of CS+ trials and punishment on 50% of CS-
% trials
rando = rand(n,1);
r = zeros(n,1);
r(cueType == 1 & ~reversed & rando > 0.2) = rewardSize;
r(cueType == 1 & reversed & rando > 0.5) = punishSize;
r(cueType == 2 & reversed & rando > 0.2) = rewardSize;
r(cueType == 2 & ~reversed & rando > 0.5) = punishSize;

% run model
model = kalmanRW(X,r,param);
% collect data
rhat = zeros(n,1); % Predicted reward
pe = zeros(n,1); % prediction error
w = zeros(2,n); % weights
Kn = zeros(2,n); % Kalman gain
offDiag = zeros(n,1); % off diagonal term in posterior weight covariance matrix
onDiag = zeros(2,n); % on diagonal terms
for counter = 1:n
   rhat(counter) = model(counter).rhat;
   pe(counter) = model(counter).dt;
   w(:,counter) = model(counter).w0;
   Kn(:,counter) = model(counter).K;
   offDiag(counter) = model(counter).C(2,1); % covariance matrix is symmetric so bottom left or top right corner of 2,2 matrix are equivalent
   onDiag(:,counter) = [model(counter).C(1,1); model(counter).C(2,2)];
end

ensureFigure('Reversal_KalmanRW', 1);
% figure;
subplot(2,2,1);
plot(trials, w');
ylim = get(gca, 'YLim');
revTrials = 0:revFreq:n;
line(repmat(revTrials, 2, 1), repmat(ylim', 1, length(revTrials)), 'Color', 'g');
ylabel('Weights (mean)');
xlabel('Trials');
subplot(2,2,2);
plot(trials, onDiag');
ylim = get(gca, 'YLim');
line(repmat(revTrials, 2, 1), repmat(ylim', 1, length(revTrials)), 'Color', 'g');
ylabel('Weights (variance)');
xlabel('Trials');
subplot(2,2,3);
plot(trials, offDiag);
ylim = get(gca, 'YLim');
line(repmat(revTrials, 2, 1), repmat(ylim', 1, length(revTrials)), 'Color', 'g');
ylabel('Weights (covariance)');
xlabel('Trials');
subplot(2,2,4);
plot(trials, pe);
ylim = get(gca, 'YLim');
line(repmat(revTrials, 2, 1), repmat(ylim', 1, length(revTrials)), 'Color', 'g');
ylabel('Prediction error');
xlabel('Trials');

saveas(gcf, fullfile(savepath, 'Reversals_KalmanRW'), 'fig');
saveas(gcf, fullfile(savepath, 'Reversals_KalmanRW'), 'jpeg');

figure; 
plot(trials,Kn);

%%
%% experiment 2-  on each trial 2 stimuli are delivered together as well as reward in first phase, then cue 1 is withheld
% shows that compound conditioning generates negative covariance between
% weight estimates owing to their linear combination to predict reward value on each trial (reward is scalar)
% I think that's the right explanation...    Saying that weight
% configuration to sum to 1 (1 = reward value) is not really correct
% because you can not deliver reward at all and it makes no difference to
% the covariance update.  I think it just has to do with the fact that
% reward is a scalar and Reward = H * Weights + measurement covariance;
% the transformation matrix H that maps weights onto reward corresponds to
% [1; 1] initially meaning that both cue 1 and cue 2 are present....
param = KTD_defparam;
param.s = 1; 

n = 40; % 10 trials
X = [ones(n/2,2)];
X = [X; [zeros(n/2,1) ones(n/2,1)]];

r = [ones(n,1)];
offdiag = zeros(n,1);
ondiag = zeros(n,2);
model = kalmanRW(X,r, param);
PE = zeros(n,1);
weights = zeros(n,2);
for counter = 1:n
    offdiag(counter) = model(counter).C(2,1);
    ondiag(counter, :) = diag(model(counter).C)';
    PE(counter) = model(counter).dt;
    weights(counter,:) = model(counter).w0'; 
    Kn(:,counter) = model(counter).K;    
end
% throw in prior covariance
Cprior = param.c*eye(2) + param.q*eye(2);
ensureFigure('Compound_conditioning', 1);
% figure; 
subplot(2,2,1);
plot(0:n, [0; offdiag], 'g'); hold on;
plot(0:n, [[Cprior(1,1) Cprior(2,2)]; ondiag]);
ylabel('Weight covariance');
legend({'off diag', 'on diag, cue 1', 'on diag, cue 2'});
xlabel('trial #');
subplot(2,2,2);
plot(PE);
ylabel('prediction error');
xlabel('Trials');
subplot(2,2,3);
plot(weights);
ylabel('weights');
xlabel('Trials');

saveas(gcf, fullfile(savepath, 'Compound_KalmanRW'), 'fig');
saveas(gcf, fullfile(savepath, 'Compound_KalmanRW'), 'jpeg');
