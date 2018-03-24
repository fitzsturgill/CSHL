

% math is weird
x = linspace(0,20,100);

y = exp(1i * x);
y2 = exp(x);
ensureFigure('test', 1);
subplot(2,2,1); plot(y, '-'); title('unit circle');
subplot(2,2,2); plot(x, real(y), '-'); title('real');
subplot(2,2,3); plot(x, imag(y), '-'); title('imaginary');
subplot(2,2,4); plot(x, angle(y), '-'); title('angle');