%% experiment-  on each trial 2 stimuli are delivered together as well as reward
param = KTD_defparam;
param.s = 1; % make reward variance large so that model learns slowly

n = 1000; % 10 trials
X = [ones(n,2)];

r = [ones(n,1)];
offdiag = zeros(n,1);
ondiag = zeros(n,1);
model = kalmanRW(X,r, param);
PE = zeros(n,1);
for counter = 1:n
    offdiag(counter) = sum(model(counter).C([2;3]));
    ondiag(counter) = sum(diag(model(counter).C));
    PE(counter) = model(counter).dt;
end
figure;
plot(1:n, offdiag, 'r'); hold on;
plot(1:n, ondiag, 'g');
ylabel('Weight covariance');
legend({'off diag', 'on diag'});
xlabel('trial #');

figure;
plot(PE);


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
