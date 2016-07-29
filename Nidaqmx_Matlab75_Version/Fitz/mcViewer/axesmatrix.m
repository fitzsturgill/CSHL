function hax = axesmatrix(nrow,ncol,ind,params,fid)
% function hax = axesmatrix(nrow,ncol,ind,params,fid)
%
% function is similar to subplot but allow more control of the
% 1) postion the matrix of axes within the figure
% 2) spacing between individual axes within the axes (or subplot) matrix
%
% INPUT
%              - nrow, number of rows in axis matrix
%              - ncol, number of columns in axis matrix
%              - ind, index of location in axis matrix (starting from
%              top-left) (can be a vector)
%              - params (optional) struct defines the "matrix of subplots"
%                  params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].
%                       LEFT and TOP are offset from top left corner of figure.
%                       WIDTH and HEIGHT define size of axesmatrix as fraction of the figure (after matmargin and figmargin are subtracted see below)
%                       e.g. [0.25 0 0.5 0.5] will result in a plot matrix starting a 1/4 of the way down the figure on the far left size and occupy a width and height of half the figure
%                  params.figmargin (optional) margins within axesmatrix: [ LEFT RIGHT TOP BOTTOM] % default [0 0 0 0]
%                  params.matmargin (optional) margins within whole figure: [ LEFT RIGHT TOP BOTTOM] % default [0 0 0 0]
%                  params.cellmargin (optional) margin within each cell (defines offset and size of each axes with a cell) % default [0 0 0 0]
%              - fid (optional), figure handle where hax will be created
%
% USAGE:
%  params.matpos = [0 0 0.5 0.75]
%  row = 2; col 5;
%  hax =axesmatrix(params,row,col)
%
% BA 011810
%

% *TO DO
% add do not create new axes if identical axes exists.
% add support for plots that span multiple cells
def_matpos = [0 0 1 1];
def_figmargin = [0.03 0.03  0.025  0 ];  % [ LEFT RIGHT TOP BOTTOM]
def_matmargin = [0 0  0  0 ];
def_cellmargin = [0.01 0.01  0.04  0.04 ];
if nargin < 5 || isempty(fid);fid = gcf;   end

if nargin < 4 || isempty(params); params = struct([]) ;   end

if ~isfield(params,'matpos');    params(1).matpos = def_matpos;end
if ~isfield(params,'figmargin');    params(1).figmargin = def_figmargin;end
if ~isfield(params,'matmargin');    params(1).matmargin = def_matmargin;end
if ~isfield(params,'cellmargin');    params(1).cellmargin = def_cellmargin;end

%         calculate [width height] of each cell in matrix
MATE = [(params.matpos(3)*1-(sum(params.figmargin ([1:2]))+sum(params.matmargin ([1:2]))))/ncol...
    (params.matpos(4)*1-(sum(params.figmargin ([3:4]))+sum(params.matmargin ([3:4]))))/nrow];

for i = 1:length(ind)
    % calculate C and R from index
    R = ceil(ind(i)/ncol);
    C = mod(ind(i),ncol)+ ~mod(ind(i),ncol)*ncol;
    
    if R > nrow || C> ncol; error('ind exceeds matrix size'); end
    
    axpos = [ (C-1)*MATE(1)+params.matpos(1)+params.cellmargin(1)+params.figmargin(1)+params.matmargin(1)...
        1-(params.figmargin(3)+params.matpos(2)+params.matmargin(3)+(R)*MATE(2)-params.cellmargin(4)) ...
        MATE(1)-params.cellmargin(2) MATE(2)-params.cellmargin(3)];
    
    hax(i) = axes('Position',axpos,'Parent',fid);
end