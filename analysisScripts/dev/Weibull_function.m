

%% exploring weibull distribution, CDF form

x = 0:100;
alpha = [0.5 1:5:26];
beta = [0.5 1:5:26];
% weibull = 1 - e^(-1 * (x/alpha)^beta)
ensureFigure('weibull', 1);
subplot(1,2,1); hold on; title('vary beta');
subplot(1,2,2); hold on; title('vary alpha');
for counter = 1:length(alpha)
    % vary beta
    subplot(1,2,1);
    a = 50;
    b = beta(counter);    
    betas{counter} = num2str(b);
    plot(x, 1 - exp(-1 * (x ./ a).^b));
    % vary alpha
    subplot(1,2,2);
    a = alpha(counter);
    b = 4;    
    alphas{counter} = num2str(a);
    plot(x, 1 - exp(-1 * (x ./ a).^b));
end

subplot(1,2,1); legend(betas); legend boxoff;
subplot(1,2,2); legend(alphas); legend boxoff;

