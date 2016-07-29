function imageOut = collapse(Image, directionOfCollapse, type)
global state gh

% This function will tahe an image as array (3d stack of intensity images) and 
% will take a maximum projection along the directionOfCollapse input (either XY, XZ, or YZ);
% type can be 'average' or none

if nargin < 3
    type = '';
end
cl=class(Image);
size = size(Image);

Rows = size(1,1);
Columns = size(1,2);

if ndims(Image) == 3
	NumberOfFrames = size(1,3);
else
	NumberOfFrames=1;
	imageOut = Image;
	return
end

switch directionOfCollapse
case 'XY' % Collapse Along Z
    if strcmp(type,'average') 
        Image = mean(Image,3);
    else
        Image = max(Image,[],3);
    end
    imageOut = Image;	
case 'XZ' % Collapse along Y
     if strcmp(type,'average') 
        Image = mean(Image,1);
    else
        Image = max(Image,[],1);
    end
	imageOut = reshape(Image, Columns, NumberOfFrames);
 	imageOut = imageOut';	
case 'YZ' % Collapse along X
	 if strcmp(type,'average') 
        Image = mean(Image,2);
    else
        Image = max(Image,[],2);
    end
	imageOut = reshape(Image, Rows, NumberOfFrames);
otherwise
	disp('The direction of projection must be ''XY'' , ''XZ'', or ''YZ''.');
end
eval(['imageOut=' cl '(imageOut);']);

