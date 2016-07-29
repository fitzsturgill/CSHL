
function out1 = CFDread2(filename, num, chan)
% CFDread2 reads image files made by confoc4 or cfcNT1.1
%	 IM = CFDread(...) returns the image or stack of images
%	 Note, as of 1/18 this is a pilot program that disentangles
%	dual channel confoc data
%  first image, num=1; first channel, chan=1
%
%   Class Support
%   -------------
%   The out image is of class uint8.
%		
%   Karel Svoboda 1/10/00 Matlab 5.3
%	 svoboda@cshl.org
%

global head_info extension;

switch nargin
case 0
     disp(['-------------------------------']);  
     disp([' CFDread2 v',version])
     disp(['-------------------------------']);
     disp([' usage: CFDread2(sweep, num, chan)'])
     disp(['        ']);
	  error(['### need filename to proceed']); 
case 1
   num=1;
   chan=1;
case 2
   chan=1;
case 3
otherwise
     disp(['-------------------------------']);  
     disp([' CFDread2 v',version])
     disp(['-------------------------------']);
     disp([' usage: CFDread2(sweep, num, chan)'])
     disp(['        ']);
	  error(['### too many parameters to proceed']); 
end

if (num<1 | num > head_info.n_images); error(['### CfdRead2: ', filename, ' im num out of range.']); end
if (chan>1 & head_info.chans==1); error(['### CfdRead2: ', filename, ' only ingle channel data.']); end

%filename=sprintf('%s%03d.cfd',strcat(directory,base),sweep)

%head_info=cf4header(filename);

if head_info.chans == 3 %i.e. dual channel acquisition
   f_length=2*head_info.pixels_xy(1)*head_info.pixels_xy(2);
   xx = freadu8(filename,head_info.hdrsize+f_length*(num-1), f_length);
   xx = reshape(xx, 2, f_length/2);
   if chan==1
      xx = reshape(xx(1,:), head_info.pixels_xy(1), head_info.pixels_xy(2));
   else
      xx = reshape(xx(2,:), head_info.pixels_xy(1), head_info.pixels_xy(2)); 
   end
else  %i.e. single channel acquisition
   f_length=head_info.pixels_xy(1)*head_info.pixels_xy(2);
   xx = freadu8(filename,head_info.hdrsize+f_length*(num-1), f_length);
   xx = reshape(xx, head_info.pixels_xy(1), head_info.pixels_xy(2));
end

out1=xx;