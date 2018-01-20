%% multiply 2 gaussians 

gf = @(x,m,v) 1/sqrt(2*pi*v) * exp(-(x-m).^2 / (2 * v));
x = linspace(0,14,10000);
gf1 = gf(x, 5, 4);
gf2 = gf(x, 5, 4);
figure; plot(x, gf1, 'r'); hold on; plot(x, gf2, 'g--'); plot(x, gf1 .* gf2 * 5, 'b');
% 
% subplot(1,2,2);
% gfuse = @(x,m,v,m2,v2) (1/sqrt(2*pi*v) * exp(-(x-m).^2 / (2 * v))) .* (1/sqrt(2*pi*v2) * exp(-(x-m2).^2 / (2 * v2)));
% plot(x, gfuse(x,5,4,10,16));