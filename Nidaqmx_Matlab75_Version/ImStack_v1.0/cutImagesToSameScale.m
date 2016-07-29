function out = cutImagesToSameScale(img1,img2)
global r c
% This function will clip the 2 images if they are not the same size.
% The second one must be larger than the first.
% The larger one will be clipped, and returned as output.

if nargin~=2
	error('cutImagesToSameScale: must input 2 images.');
end

s1=size(img1);
s2=size(img2);
[r,c]=find(img2==0);
if length(s1) ~= length(s2) | length(s1) > 2 | length(s2) > 2
	error('cutImagesToSameScale: images must be 2D and the same size.');
end

if all(s1==s2)
	out=img2;
	disp('cutImagesToSameScale:Images were the same size.');
	return
end

diff=floor(.5*(s2-s1));   % Compute leftover pixels...
if any(diff<0)
	error('cutImagesToSameScale: base image smaller than registered. Select more points to align.');
end

out=[];
out=img2(diff(1)+1:(end-diff(1)),diff(2)+1:(end-diff(2)));

sO1d=size(img1);
sOut=size(out);

if all(sOut == sO1d)% Check how we did...
	tag='';
	tag=inputname(2);
	disp(['cutImagesToSameScale: Cut pixels from image ' tag '.']);
	return
else
	s1=size(img1);
	s2=size(out);
	diff=round(.5*(s2-s1));   % Compute leftover pixels...
	out=out(diff(1):(end-diff(1)),diff(2):(end-diff(2)));
end



