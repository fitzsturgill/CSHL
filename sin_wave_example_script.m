


%%
t = 0:1/100:2 - 1/100;
freq = 8;
y = sin(2 * pi * freq * t)*10;
y90 = sin(2 * pi * freq * t + pi)*10;

ensureFigure('test', 1);

plot(t, y, 'b'); hold on;
plot(t, y90, 'r');