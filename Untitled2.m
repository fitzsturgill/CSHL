% why are the slopes different?  WTF?

x = linspace(1,10,10000);
y = x;
x = x(:); y = y(:);
y2 = x;
nf = 4;
% y2 = y + randn(size(x)) * nf;
x2 = x + randn(size(x)) * nf; 
m = 5;
y = y * m;
y2 = y2 * 5;
 y2 = y2 + randn(size(x)) * nf;


figure; 
ax = [];
ax(1) = subplot(1,2,1);
scatter(x,y, '.'); hold on;
fo = fit(x,y, 'poly1');
plot(fo, 'predfunc');
legend('off');
coeffs = coeffvalues(fo);
rval = corr(x,y);
textBox(sprintf('R=%.2f, y=%.2fx + %.2f', rval, coeffs(1), coeffs(2)));%['R=' num2str(corr(x,y))]);
ax(2) = subplot(1,2,2);
scatter(x2,y2, '.'); hold on;
fo = fit(x2,y2, 'poly1');
plot(fo, 'predfunc');
legend('off');
coeffs = coeffvalues(fo);
rval = corr(x2,y2);
% textBox(['R=' num2str(corr(x2,y2))]);
textBox(sprintf('R=%.2f, y=%.2fx + %.2f', rval, coeffs(1), coeffs(2)));%['R=' num2str(corr(x,y))]);
sameXScale(ax);
sameYScale(ax);