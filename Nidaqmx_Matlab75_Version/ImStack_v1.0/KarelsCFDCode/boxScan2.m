function out=boxScan2(filename)
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

[x,y]=ginput(2);
x=(sort(x));
y=(sort(y));
rectangle('Position',[x(1),y(1),x(2)-x(1),y(2)-y(1)], 'EdgeColor','y')
%si=size(im)

timeseries(1:head_info.n_images)=0;
x=round(x);
y=round(y);

for i=1:head_info.n_images
   xx=CFDread2(filename, i);
   xx=xx(x(1):x(2),y(1):y(2));
   timeseries(i)=mean(mean(xx));
end
%pause;
out=timeseries;