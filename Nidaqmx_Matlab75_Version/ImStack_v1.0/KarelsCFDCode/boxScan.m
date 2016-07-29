function out=boxScan(filename)
% boxScan computes the average pix value in a user-define rectangle 
% inputs: 	
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
     disp([' boxScan v',version])
     disp(['-------------------------------']);
     disp([' usage: boxScan(filename)'])
     error(['### not enough parameters to proceed']); 
end

if head_info.n_images < 2
   disp(['-------------------------------']);  
   disp([' boxScan v',version])
   disp(['-------------------------------']);
   disp([' usage: boxScan(filename)'])
   error(['### file contains only a single image']); 
end

[x1,y1]=ginput(1)
l1=line([x1 head_info.pixels_xy(1)], [y1 y1], 'Color','y');
l2=line([x1 x1], [y1 head_info.pixels_xy(2)], 'Color','y');
%x=(sort(x));
%y=(sort(y));
[x2,y2]=ginput(1)
set(l1, 'LineStyle', 'none');
set(l2, 'LineStyle', 'none');
rectangle('Position',[x1,y1,x2-x1,y2-y1], 'EdgeColor','y')
%si=size(im)

timeseries(1:head_info.n_images)=0;
x1=round(x1);
y1=round(y1);
x2=round(x2);
y2=round(y2);

for i=1:head_info.n_images
   xx=CFDread2(filename, i);
   xx=xx(x1:x2,y1:y2);
   timeseries(i)=mean(mean(xx));
end
%pause;
out=timeseries;