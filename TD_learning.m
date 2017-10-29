%%
ensureFigure('TD_learning', 1); axes; hold on;
nPoints = 20;
v = zeros(1,nPoints); % discounted value
r = zeros(1,nPoints); % reward
rTime = round(0.8 * nPoints);
rP = 0.5; % reward probability
x = zeros(1,nPoints); % stimulus
x(round(nPoints * 0.2):end) = 1; % stimulus on at point 20
w = zeros(1,nPoints - 1); % weights
pe = zeros(1,nPoints - 1); % prediction error td
alpha = 0.2; % learning rate
gamma = 0.95; % discount factor 0 < gamma < 1
trials = 1000;
cmap = parula;
for i = 1:trials
    thisReward = r;
    if rand(1) >= rP
        thisReward(rTime) = 1;
    end
    % pe = r(t) + gamma * v(t+1) - v(t)
    pe = thisReward(1:end-1) + gamma * v(2:end) - v(1:end-1);
    w = w + alpha * x(1:end-1) .* pe; % adjustment of weights
    v(1:end-1) = x(1:end-1) .* w; % update of discounted value function
    if rem(i,10) == 0
        plot(v, 'Color', cmap(round(i/trials * size(cmap, 1)), :));
    end
end
plot(v);
stem([0 diff(x)],'Color', 'black', 'Marker', 'none')
plot(r,'r')