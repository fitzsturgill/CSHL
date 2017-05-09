%%
 rawTimes = sessions(1).SessionData.RawEvents.Trial{1}.Events.Port3In;
 wheelTimes = rawTimes(rawTimes < 30);
 wheelY = ones(size(wheelTimes));
 wheelTimes = [0 wheelTimes 30];
 wheelY = cumsum([0 wheelY 0]);
 
 [wheelY_new, wheelTimes_new] = resample(wheelY, wheelTimes, 20, 'linear');
 
 ensureFigure('test', 1);
 plot(wheelTimes, wheelY); hold on;
 plot(wheelTimes_new, wheelY_new);
 
 
 %%
D = 12.7; % diameter of wheel in cm
pr = 200; % pulses/rotation 
dpp = pi * 12.7 / 200;
rawTimes = sessions(1).SessionData.RawEvents.Trial{5}.Events.Port3In;
wheelTimes = rawTimes(rawTimes < 30);
edges = 0:1/20:30;
% wheel_X = cumsum(histcounts(rawTimes, edges) * dpp);
wheel_X = smooth(cumsum(histcounts(rawTimes, edges)), 20);  
wheel_t = 0:1/20:30-1/20;
  
ensureFigure('test', 1);
subplot(2,1,1);
plot(wheelTimes, cumsum(ones(size(wheelTimes)) * dpp), 'b'); hold on;
plot(wheel_t, wheel_X, 'g');

subplot(2,1,2);
plot(wheel_t(1:end-1), diff(wheel_X), 'b'); hold on
plot(wheel_t, gradient(wheel_X), 'g');
 