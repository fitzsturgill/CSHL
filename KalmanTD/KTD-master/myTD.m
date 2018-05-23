
output = [];
%starting variables part1
phases = 2;
trials = 10;
timeSteps = 10;
totalTS = trials * timeSteps * phases;
numCues = 2;

%establish variables part2
disc = 0.98; % temporal discount factor
diffVar = .01; % diffusion variance (chaos)
generalDistrust = 1; %reduces kalman gain by how much we distrust our data (noise)
initialCov = 1;
perTrialCues = zeros(timeSteps,numCues);

%phase1
perTrialCues(3:6,1) = 1;
cueData = [];
for i = 1:numCues
    cueData = [cueData diag(perTrialCues(:,i))];
end
cueData = repmat(cueData,trials,1);

perTrialReward = zeros(timeSteps,1);
perTrialReward(6) = 1;
rewardData = repmat(perTrialReward,trials,1);

%phase2
perTrialCues = zeros(timeSteps,numCues);
perTrialCues(3:6,2) = 1;
cueData2 = [];
for i = 1:numCues
    cueData2 = [cueData2 diag(perTrialCues(:,i))];
end
cueData2 = repmat(cueData,trials,1);

perTrialReward = zeros(timeSteps,1);
perTrialReward(6) = 1;
rewardData2 = repmat(perTrialReward,trials,1);


rewardData = [rewardData; rewardData2;];
cueData = [cueData; cueData2;];
cueData = [cueData; zeros(1,numCues * timeSteps)];
%establish time chains

weightData = zeros(numCues * timeSteps,1);


%final setup
priordCov = diffVar * eye(numCues * timeSteps);
Cov = initialCov * eye(numCues * timeSteps);
K = zeros(totalTS,numCues);


kalman

%run TD
for i = 1:totalTS
    dT = cueData(i,:) - disc * cueData(i+1,:); %the temporal difference
    estReward = dT * weightData; %estimated reward
    perr = rewardData(i) - estReward; %self-explanatory prediction error calculation
    Cov = Cov + priordCov; %prior covariance update
    K = (Cov * dT') / ((dT * Cov * dT') + generalDistrust);
    w0 = weightData;
    weightData = weightData +  K*perr;
    Cov = Cov -  K*dT*Cov; %posterior covariance update
    
    output(i).w = weightData;
    output(i).w0 = w0;
    output(i).C = Cov;
    output(i).K = K;
    output(i).dt = perr;
    output(i).rhat = estReward;
end
 

figure;
%plot(output();

