function out=lineScan(filename)
% LINESCAN computes the average pix value between user-defined rows
%		
%   
%   Class Support
%   -------------
%   
%		
%   Karel Svoboda 1/13/00 Matlab 5.3
%	 svoboda@cshl.org
%
global head_info

if nargin < 1 
     disp(['-------------------------------']);  
     disp([' MOVIESCAN v',version])
     disp(['-------------------------------']);
     disp([' usage: MOVIESCAN(im, [x1,x2],[y1,y2])'])
     error(['### not enough parameters to proceed']); 
end

[x1,y1]=ginput(1);
%whos x1
line([1 head_info.pixels_xy(1)], [y1 y1], 'Color','y');
%pause
[x2,y2]=ginput(1);
line([1 head_info.pixels_xy(1)], [y2 y2], 'Color','y');

%x=(sort(x));
%y=(sort(y));
%rectangle('Position',[x(1),y(1),x(2)-x(1),y(2)-y(1)], 'EdgeColor','y')
%line([1 32], [y(2) y(2)], 'Color','y');

timeseries(1:(head_info.n_images*head_info.pixels_xy(1)))=0;

y1=round(y1);
y2=round(y2);
for i=1:head_info.n_images
   xx=CFDread2(filename, i);
   xx=xx(y1:y2,:);
   timeseries(((i-1)*head_info.pixels_xy(1)+1):((i)*head_info.pixels_xy(1)))=mean(xx);
end
%pause;
out=timeseries;