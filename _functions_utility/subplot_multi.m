function [figure_handle,subplot_handles] = subplot_multi(...
    x_in, ...               %1
    y_in,...                %2
    plot_varargin,...       %3
    x_y_labels,...          %4
    plot_legends,...        %5
    subplot_shape,...       %6
    subplot_margins,...     %7
    subplot_varargin)       %8


%   subplot_multi() will use two intputs (x_in,y_in) to build a series of
%   plots corresponding to each row of x_in,y_in, in a variety of subplots.
%   By default, subplots are generated automatically in square(ish) format.
%
%   INPUTS
%
%   Note: to skip any input variable, provide a dummy input with {} or [].
%
%       x_in and y_in
%           Data to pbe plotted. Should be of the format size(x_in) =[num_lines,line_length]
%           That is, rows of x_in are individual lines, and columns are points in
%           each line.
%
%       plot_varargin
%           If you pass a cell array {}, the entries in it will be passed
%           to plot() to alter the figure composition. For example, if you
%           want all plots to have thick blue lines, pass
%           {'b','LineWidth',2}. Inputs here follow the convention of
%           Matlab's plot() function.
%
%       legend
%           A cell array of strings for legends.  Note that this is
%           expanded to allow for only 1 string per subplot! (Since that is
%           all we can plot anyway.)
%
%       m,n
%           By default we will attempt a square-ish number of subplots.  TO
%           over-ride this, specify subplot rows (m) and columns (n)
%
%       subplot_margins
%           Following the convention of subplot_tight, define the vertical
%           and horizontal margins between subplots as
%           [vertical,horizontal].
%           NOTE :: USE NaN to force use of normal subplot() function
%
%       subplot_varargin
%           If you pass a cell array {}, the entries in it will be passed
%           to subplot() to alter subplot parameters.  See subplot() for
%           details.
%
%   OUTPUTS
%       
%       subplot_handles :: an array of subplot handles
%
% Alex Vaughan, 2015

% Data is helpful.
if nargin < 2,
    error('Must provide at least x_in and y_in')
end
n_subplots = size(y_in,1);
if size(x_in,1) == 1,
    x_in = repmat(x_in,n_subplots,1);
end
size(x_in)
size(y_in)
assert(all(size(x_in)==size(y_in)),'x_in and y_in do not appear to be compatible sizes.')

% Plot varargin
if nargin < 3 || numel(plot_varargin) == 0
    plot_varargin = {};
end

if nargin < 4 || numel(x_y_labels) == 0
   x_labels = {};
   y_labels = {};
else
    x_labels = x_y_labels{1};
    y_labels = x_y_labels{2};
end

% Cell array of strings for the legend.
if nargin < 5 || numel(plot_legends) == 0
    plot_legends = {};
end
% Expand to fill a full cell array.
assert(numel(plot_legends) <= 1 || numel(plot_legends) == n_subplots,'The length of cell array plot_legends must be either 0, 1, or equal to the number of plots')
if length(plot_legends) == n_subplots
    % One string per subplot
    all_legends = plot_legends;
else
    % One string for ALL subplots
    [all_legends{1:n_subplots}] = deal(plot_legends);
end

% Define row and column sizes for the subplots
if nargin < 6 || numel(subplot_shape) == 0
    subplot_rows = ceil(sqrt(n_subplots));
    subplot_columns = ceil(n_subplots/subplot_rows);
else
    subplot_rows = subplot_shape(1);
    subplot_columns = subplot_shape(2);
end

% Margin width.
if nargin >=7 && isnumeric(subplot_margins) && isnan(subplot_margins)
    use_tight_subplot = 0;
else
    use_tight_subplot = 1;
    if nargin < 7 || numel(subplot_margins) == 0,
        subplot_margins = [0.05,0.05];
    end
end

% Elements of subplot_varargin are passed to subplot() as varargin inputs.
if nargin < 8,
    subplot_varargin = {};
end

%% Perform actual plotting
figure_handle = figure;  clf
subplot_handles = zeros(n_subplots,1);
for this_subplot = 1:n_subplots,
    
    % Generate subplot
    if use_tight_subplot
        subplot_handles(this_subplot) = subplot_tight(subplot_rows,subplot_columns,this_subplot,subplot_margins,subplot_varargin{:});
    else
        subplot_handles(this_subplot) = subplot(subplot_rows,subplot_columns,this_subplot,subplot_varargin{:});
    end
    
    % Perform the actual plotting.
    plot(x_in(this_subplot,:),y_in(this_subplot,:),plot_varargin{:});
    
    % Legends and xlabel/ylabel.  These are mostly nested cells.
    if ~isempty(plot_legends)
        legend(all_legends{this_subplot});
    end
    if length(x_labels) >= n_subplots
        xlabel(x_labels{this_subplot});
    end
    if length(y_labels) >= n_subplots
        ylabel(y_labels{this_subplot});
    end
    
    
end


