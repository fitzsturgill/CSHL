                
%%        commented code below is snippet to show what different coefficients do to an exponential
% note! c in the example below is the inverse of the time constant, a is
% the asymptote at x = Inf, a + b is the y intercept
x = 1:1000;

a = 100;
b = 1;
c = -.01;

ensureFigure('test', 1);
subplot(2,2,1);
plot(x, a + b*exp(c * x), 'k'); hold on;
plot(x, a + b * 2 *exp(c * x), 'b');
plot(x, a + b * exp(c * 2 * x), 'r');
plot(x, a * 1.01 + b * exp(c * x), 'y');

%
c = 0.01;
subplot(2,2,2);
plot(x, a + b*exp(c * x), 'k'); hold on;
plot(x, a + b * 2 *exp(c * x), 'b');
plot(x, a + b * exp(c * 2 * x), 'r');
plot(x, a * 1.01 + b * exp(c * x), 'y');
title('c positive');
%
c = -0.01;
b = -1;
subplot(2,2,3);
plot(x, a + b*exp(c * x), 'k'); hold on;
plot(x, a + b * 2 *exp(c * x), 'b');
plot(x, a + b * exp(c * 2 * x), 'r');
plot(x, a * 1.01 + b * exp(c * x), 'y');
title('b negative');