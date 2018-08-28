function plot_row_combinations(...
    data_in,        ... % REQUIRED
    plot_depth,     ... % optional
    plot_varargin,...   % optional
    x_y_labels,...      % optional
    plot_legends,...    % optional  
    subplot_shape,...   % optional  
    subplot_margins,... % optional  
    subplot_varargin)   % optional  

%
%   plot_row_combinations will plot ~all pairwise combinations of the rows
%       of data_in.  This is useful, for example, for plotting individual
%       principal components against each other.
%
%   INPUTS
%       
%       data_in :: size(n_dimensions, n_samples).  By default, we will plot
%           n^2 plots, with each dimension contrasted with each other
%           dimension.  
%
%       plot_depth :: Limits plots to the pairwise combinations of the
%           first plot_depth rows of data_in.  Total number of plots is
%           then (plot_depth^2).
%
%       NOTE :: 
%       DEFAULTS ::  By default, will plot random data.
%       
%   OUTPUTS
%       
%       none.
%
% Alex Vaughan, 2015
%

% Default data
if nargin < 2,
    warning('Using default data for plot_row_combinations!')
    data_in = rand(5,200);
    plot_depth = 5;
end

% DEFAULTS ARE PASSED TO SUBPLOT_MULTI
if nargin < 3
    plot_varargin = {'k.'};
end
if nargin < 4
    x_y_labels = {};
end
if nargin < 5
    plot_legends = {};
end
if nargin < 6
    subplot_shape = {};
end
if nargin < 7
    subplot_margins = {};
end
if nargin < 8
    subplot_varargin = {};
end

% Setup data and labels
n_dimensions = min(size(data_in,1),plot_depth);

% data_to_plot provides indices into the rows of data_in to allow plotting 
% each row against each other 
[data_to_plot{1:2}] = ind2sub([n_dimensions,n_dimensions],1:(n_dimensions^2));

% xlabel and ylabel are given labels based on the indices in data_to_plot
x_y_labels{1} = arrayfun(@(x) {sprintf('Dimension %.0f',x)},data_to_plot{1},'UniformOutput',false);
x_y_labels{2} = arrayfun(@(y) {sprintf('Dimension %.0f',y)},data_to_plot{2},'UniformOutput',false);

%% Call subplot_multi to do the plotting.

[figure_handle,subplot_handles] = subplot_multi(...
    data_in(data_to_plot{1},:),...      % 1 :: x_in
    data_in(data_to_plot{2},:),...      % 2 :: x_out
    plot_varargin,...                   % 3 :: plot_varargin
    x_y_labels,...                      % 4 :: x_y_labels
    plot_legends,...                    % 5 :: plot_legends
    subplot_shape,...                   % 6 :: subplot_shape
    subplot_margins,...                 % 7 :: subplot_margins - use NaN to force use of normal subplot function
    subplot_varargin);                  % 8 :: subplot_varargin



