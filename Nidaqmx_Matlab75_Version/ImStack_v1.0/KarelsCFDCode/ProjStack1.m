
function out=ProjStack1(filename, chan, range, dim, theta, padding)
% ProjStack reads a 3D image stack and computes arbitrary Max projections
%				
%   Karel Svoboda 3/10/00 Matlab 5.3
%	 svoboda@cshl.org
%
% INPUT PARAMETERS
% filename: full file name i.e. 'd:\data\k001a001.cfd'
% chan: optional (default 1), range 1,2
% range: optional (default [1,max]), two-element vector where [range(1) range(2)]=[min max]
% dim: optional (default 'z'); range x', y', z'
% 
% OUTPUT PARAMETERS
% 2-dimensional uint8 array

global head_info;

switch nargin
case 0
     disp(['-------------------------------']);  
     disp([' ProjZ v',version])
     disp(['-------------------------------']);
     disp([' usage: ProjZ(filename, chan, range)'])
     disp(['        ']);
	  error(['### need filename to proceed']); 
case 1
   range=[1, head_info.n_images];
   chan=1;
   dim='z';
   theta=0;
case 2
   range=[1, head_info.n_images];
   dim='z';
	theta=0;
case 3 
   dim='z'; 
   theta=0;
case 4   
   theta=0
case 5
otherwise
     disp(['-------------------------------']);  
     disp([' ProjZ v',version])
     disp(['-------------------------------']);
     disp([' usage: ProjZ(filename, chan, range)'])
     disp(['        ']);
	  error(['### too many parameters to proceed']); 
end

range
%force range to make sense i.e. ramge(2)>range(1); range (1) > 0; range(2)<n_images
if range(1) > range(2); range=range([2,1]);end
if range(1) < 1; range(1)=1; end
if range(2) > head_info.n_images; range(2)=head_info.n_images; end

if (chan>1 & head_info.chans==1); error(['### CfdRead2: ', filename, ' only single channel data.']); end

if head_info.chans == 3 %i.e. dual channel acquisition
   f_length=2*head_info.pixels_xy(1)*head_info.pixels_xy(2)*(range(2)-range(1)+1);
   offset=head_info.hdrsize+2*head_info.pixels_xy(1)*head_info.pixels_xy(2)*(range(1)-1);
   xx = freadu8(filename,offset, f_length);
   xx = reshape(xx, 2, f_length/2);
   if chan==1
		xx = reshape(xx(1,:), head_info.pixels_xy(1), head_info.pixels_xy(2), (range(2)-range(1)+1));
	else
		xx = reshape(xx(2,:), head_info.pixels_xy(1), head_info.pixels_xy(2), (range(2)-range(1)+1)); 
	end
else  %i.e. single channel acquisition
   f_length=head_info.pixels_xy(1)*head_info.pixels_xy(2)*(range(2)-range(1)+1);
   offset=head_info.hdrsize+head_info.pixels_xy(1)*head_info.pixels_xy(2)*(range(1)-1);
   xx = freadu8(filename, offset, f_length);
   xx = reshape(xx, head_info.pixels_xy(1), head_info.pixels_xy(2),(range(2)-range(1)+1));
end

switch dim
	case 'z'
      xx=permute(xx,[3 2 1]);
   case 'y'
      xx=permute(xx,[2 1 3]);
   case 'x'
      if theta ~= 0;
         temp=imrotate(xx(:,:,1),theta);
         si=size(temp);
         MaxVal(1:si(1)*(range(2)-range(1)+1))=0;
         MaxVal(1:si(1))=max(temp);
         for i=2:(range(2)-range(1)+1)
            temp=imrotate(xx(:,:,i),theta);
            MaxVal((si(1)*i+1):(si(1)*(i+1)))=max(temp);
         end
         MaxVal=reshape(MaxVal(1:si(1)*(range(2)-range(1)+1)), [si(1) range(2)-range(1)+1]);
         MaxVal=uint8(MaxVal);
      end
   otherwise
   end
   
if theta==0; MaxVal=max(xx); end
MaxVal=squeeze(MaxVal);
out=MaxVal';