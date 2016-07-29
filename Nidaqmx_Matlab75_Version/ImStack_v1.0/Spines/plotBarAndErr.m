function center = plotBarAndErr(b,errdata)
global state gh

h = bar(b);
xdata = get(h, 'XData');
ydata = get(h ,'YData');

% Determine number of bars
sizz = size(b);
nb = sizz(1)*sizz(2);
xb = [];
yb = [];
xb = [xb xdata];
yb = [yb ydata];


% To find the center of each bar - need to look at the output vectors xb, yb
% find where yb is non-zero - for each bar there is a pair of non-zero yb
% find where yb is non-zero - for each bar there is a pair of non-zero yb
%values.
% The center of these values is the middle of the bar

nz = find(yb);
for i = 1:nb,
    center(i) = (xb(nz(i*2))-xb(nz((i*2)-1)))/2 + xb(nz((i*2)-1));
end;

% To place the error bars - use the following:

hold on;
h=errorbar(center, b, errdata);
set(h(1),'linewidth',2);            % This changes the thickness of the errorbars
set(h(1),'color','r');              % This changes the color of the errorbars
set(h(2),'linestyle','none');       % This removes the connecting
