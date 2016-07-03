function lh = linecustommarker(x,y,xymark,msiz,hAxes)
% lh = linecustommarker(x,y,xymark,msiz,hAxes)
% INPUT
%      x, y values
%      xymark marker template (default vertical line)
%      msiz optional argument to scale markersize by x and y limits of axis
%      Tip: leave empty for spike raster in order to get unit length lines
%      in y.
% from Jonas Rose:
% http://www.mathworks.in/matlabcentral/newsreader/view_thread/77718
if nargin<5 % default msizy to lines that span less than one unit in Y
hAxes = gca;
end
% create a marker template
if nargin<3 || isempty(xymark)
    %      xmark=[-1 0 1 -1];	% a triangle
    %      ymark=[1 0 1 1];
    xmark=[0 0];	% a baton
    ymark=[-1 1];
else
    xmark = xymark(1,:);
    ymark =  xymark(2,:);
end
if nargin<4 || isempty(msiz)% default msizy to lines that span less than one unit in Y
    msizy = 0.440;
    msizx = 0.001;                                    % this value doesn't seem to impact anythig about the plot
else %this the original code it allows rasterline to be changed in proportion
    %     by msiz x the size of the current axis
    % this is NOT desirable for a spike raster cause the length of raster
    % lines to changed based on total number of trials
    
%     msiz=.022;	% rel marker size 
    
    % get proper geometry...
    lh1=line(x,y);
    msizx=msiz*range(get(hAxes,'xlim'));
    msizy=msiz*range(get(hAxes,'ylim'));
    delete(lh1);
end

% prepare x/y-markers
     x=x(:);
     y=y(:);
     xx=repmat(x,1,size(xmark,2))+...
        msizx*repmat(xmark,size(x));
     xx=[xx nan(size(x))];
     xx=reshape(xx.',[],1);
     yy=repmat(y,1,size(ymark,2))+...
        msizy*repmat(ymark,size(y));
     yy=[yy nan(size(y))];
     yy=reshape(yy.',[],1);
% ...and draw them
     lh=line(xx,yy,'linewidth',0.75,'Parent',hAxes);
     
    