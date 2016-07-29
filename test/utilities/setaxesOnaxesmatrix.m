function setaxesOnaxesmatrix(hAxes,nrow,ncol,ind,params,fid)
% function setaxesOnaxesmatrix(hAxes,nrow,ncol,ind,params,fid)
% like axesmatrix, but sets existing hAxes to matrix positions rather than
% creating them
% BA
def_matpos = [0 0 1 1];
def_figmargin = [0.03 0.03  0.025  0 ];  % [ LEFT RIGHT TOP BOTTOM]
def_matmargin = [0 0  0  0 ];
def_cellmargin = [0.01 0.01  0.04  0.04 ];
if nargin < 6 || isempty(fid);fid = gcf;   end
if nargin < 4 || isempty(ind); ind = [1:prod(size(hAxes))];   end

if nargin < 5 || isempty(params); params = struct([]) ;   end

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
    
    set(hAxes(i),'Position',axpos,'Parent',fid);
end