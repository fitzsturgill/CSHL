

param = KTD_defparam;
rewardSize = 1;
punishSize = -1;
revFreq = 100; % block size (normal or reversed contingencies alternate)
n = 1000;
trials = (1:n)';
cueType = randi(2,n,1); % each trial contains either stimulus 1 or 2
X = zeros(n,2);
X(cueType == 1, 1) = 1;
X(cueType == 2, 2) = 1;


reversed = floor(trials / revFreq);
reversed = rem(reversed, 2) == 1;

rando = rand(n,1);
r = zeros(n,1);
r(cueType == 1 & ~reversed & rando > 0.2) = rewardSize;
r(cueType == 1 & reversed & rando > 0.5) = punishSize;
r(cueType == 2 & reversed & rando > 0.2) = rewardSize;
r(cueType == 2 & ~reversed & rando > 0.5) = punishSize;

% run model
model = kalmanRW(X,r,param);
% collect data
rhat = zeros(n,1);
pe = zeros(n,1);
w = zeros(2,n);
Kn = zeros(2,n);
offDiag = zeros(n,1);
onDiag = zeros(2,n);
for counter = 1:n
   rhat(counter) = model(counter).rhat;
   pe(counter) = model(counter).dt;
   w(:,counter) = model(counter).w0;
   Kn(:,counter) = model(counter).K;
   offDiag(counter) = model(counter).C(2,1); % covariance matrix is symmetric so bottom left or top right corner of 2,2 matrix are equivalent
   onDiag(:,counter) = [model(counter).C(1,1); model(counter).C(2,2)];
end

ensureFigure('Reversal_KalmanRW', 1);
subplot(2,2,1);
plot(trials, w');
ylabel('weights');
subplot(2,2,2);
plot(trials, onDiag');
subplot(2,2,3);
plot(trials, offDiag);