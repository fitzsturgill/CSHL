% trying to figure out how negative covariance arises
%%
param = KTD_defparam;
param.s = 1;
% overshadowing
n = 100;
X = [ones(n,2)];

r = [ones(n,1)];

model = kalmanRW(X,r, param);
VarW = zeros(n,1);
pe = zeros(n,1);
for counter = 1:n
    VarW(counter) = sum(diag(model(counter).C)); 
    pe(counter) = model(counter).dt;
end
ensureFigure('test', 1);
subplot(2,1,1);
plot(1:n, VarW);
subplot(2,1,2);
plot(1:n, pe);
% results.W(i,:) = model(end).w0;
% results.k(i,:) = model(end).K;

%% try eliminating weight diffusion variance and reward variance (set them to 0 like I did before when I was confused)
    param.c = 0.011;                % prior variance
    param.s = 1; % 1            % noise variance
    param.q = 0.01; % 0.01              % diffusion variance
    param.g = 0.98;             % discount factor
%     param.K = 6;                % number of microstimuli
%     param.sigma = 0.08;         % width of microstimuli
%     param.decay = 0.985;        % trace decay
    param.lr = 0.3;             % learning rate for standard TD
    param.TD = 0;               % use standard TD?
    
    % overshadowing
n = 1;
X = [ones(n,2); zeros(n,2)];

r = [ones(n,1); zeros(n,1)];

model = kalmanRW(X,r, param);

%% "Kalman gain grows monotonically with stimulus uncertainty encoded in covariance matrix diagonals"
Kn = zeros(100,2);
n = 1;
X = [ones(n,2); zeros(n,2)];
r = [ones(n,1); zeros(n,1)];
c = linspace(0, 10, 100)';
for counter = 1:100
    param = KTD_defparam;
    param.c = c(counter);
    model = kalmanRW(X,r,param);
    Kn(counter,:) = model.K;
end
ensureFigure('Kn_vs_CovPrior', 1);
plot(c, Kn);

% REALIZATION- IN OVERSHADOWING YOU GET APPROXIMATELY THE SAME VALUE REWARD
% EACH TRIAL and both stimuli are presented together.  Therefore, the
% weights of each stimuli must be negatively correlated such that they add
% up to approximately 1






    