function rwb = red_white_blue(m)
% Generates a red_white_blue colormap
% With log-scale red and blue

if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

c_gradient = linspace(0,1,floor(m/2))';
half_ones = ones(floor(m/2),1);

%Ensure center point is white
if mod(m,2)
    center = [1];
else
    center = [];
end

blue = [ half_ones; center; c_gradient(end:-1:1)];
green = [  c_gradient; center; c_gradient(end:-1:1)];
red = [ c_gradient; center; half_ones];

rwb = [red green blue];

%rwb = [1-gradient(end:-1:1) 0.5*ones(size(gradient)) 1-gradient]