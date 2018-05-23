%% experiment-  on each trial 2 stimuli are delivered together as well as reward in first phase, then cue 1 is withheld
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
end
% throw in prior covariance
Cprior = param.c*eye(2) + param.q*eye(2);
ensureFigure('compound', 1);
subplot(2,2,1);
plot(0:n, [0; offdiag], 'g'); hold on;
plot(0:n, [[Cprior(1,1) Cprior(2,2)]; ondiag]);
ylabel('Weight covariance');
legend({'off diag', 'on diag, cue 1', 'on diag, cue 2'});
xlabel('trial #');
subplot(2,2,2);
plot(PE);
ylabel('predictione error');
subplot(2,2,3);
plot(weights);
ylabel('weights');



%% experiment 2-  on each trial 1 of 2 stimuli are delivered as well as reward
param = KTD_defparam;
param.s = 1; % make reward variance large so that model learns slowly

n = 100; % 10 trials
rando = rand(n,1);
X = double([rando < 0.5 rando >= 0.5]);
r = [ones(n,1)];
offdiag = zeros(n,1);
ondiag = zeros(n,1);
model = kalmanRW(X,r, param);

for counter = 1:n
    offdiag(counter) = sum(model(counter).C([2;3]));
    ondiag(counter) = sum(diag(model(counter).C));
end
figure;
plot(1:n, offdiag, 'r'); hold on;
plot(1:n, ondiag, 'g');
ylabel('Weight covariance');
legend({'off diag', 'on diag'});
xlabel('trial #');


%%
ensureFigure('test',1);
mu = [1 0];
Sigma = [1 -.5; -.5 1];
x1 = -3:.2:3; x2 = -3:.2:3;
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
F = reshape(F,length(x2),length(x1));
surf(x1,x2,F);
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
axis([-3 3 -3 3 0 .4])
xlabel('x1'); ylabel('x2'); zlabel('Probability Density');