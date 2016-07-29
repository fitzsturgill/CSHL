function loadMathImage
global gh state
% Function will compute the necessary operations on the iomages fromt he mathGUI 
% Int he image Oroicessing environments.

set(gh.mathGUI.figure1,'Pointer','watch');

% What opeartion are we doing...
operationString = get(gh.mathGUI.operation, 'String');
if iscell(operationString)
	operation=operationString{state.imageProc.mathGUI.operation};
else
	operation = operationString;
end

% Which images??

first=state.imageProc.mathGUI.fileName1;
second= state.imageProc.mathGUI.fileName2;

image1 = double(state.imageProc.cell.currentImage{first}(:,:,state.imageProc.mathGUI.startFile1:state.imageProc.mathGUI.endFile1));
image2 = double(state.imageProc.cell.currentImage{second}(:,:,state.imageProc.mathGUI.startFile2:state.imageProc.mathGUI.endFile2));

if any(size(image1) ~= size(image2))
	error('loadMathImage: Images are not the same size');
else
	switch operation
	case 'Subtract'
		if state.imageProc.mathGUI.useValue == 1
			state.imageProc.mathGUI.modImage = image1 - state.imageProc.mathGUI.mathValue;
		else
			state.imageProc.mathGUI.modImage = image1 - state.imageProc.mathGUI.weight*image2;
		end
	case 'Add'
		if state.imageProc.mathGUI.useValue == 1
			state.imageProc.mathGUI.modImage = image1 + state.imageProc.mathGUI.mathValue;
		else
			state.imageProc.mathGUI.modImage = image1 + state.imageProc.mathGUI.weight*image2;
		end
	case 'Multiply'
		if state.imageProc.mathGUI.useValue == 1
			state.imageProc.mathGUI.modImage = image1.*state.imageProc.mathGUI.mathValue;
		else
			state.imageProc.mathGUI.modImage = (image1).*(state.imageProc.mathGUI.weight*image2);
		end
	case 'Divide'
		if state.imageProc.mathGUI.useValue == 1
			if state.imageProc.mathGUI.mathValue ~= 0
				state.imageProc.mathGUI.modImage = image1./state.imageProc.mathGUI.mathValue;
			else
				disp(['Divide by zero: set new mathValue']);
				return
			end
		else
			state.imageProc.mathGUI.modImage = (image1)./(state.imageProc.mathGUI.weight*image2);
		end
	case 'Register'
		evalin('base', 'global base_points input_points');
		global base_points input_points
		state.imageProc.mathGUI.inputImage=image2;
		state.imageProc.mathGUI.baseImage=image1;
		image1=state.imageProc.mathGUI.lut*image1./max(max(image1));
		image2=state.imageProc.mathGUI.lut*image2./max(max(image2));
		h=cpselect(image2,image1);
		return
	end
end

% for frameCounter=1:size(state.imageProc.mathGUI.modImage,3)
% 	[r,c]=find((state.imageProc.mathGUI.modImage(:,:,frameCounter)<0))
% 	state.imageProc.mathGUI.modImage(r,c,frameCounter)=0;
% end
% 
if ndims(state.imageProc.mathGUI.modImage) <= 2
	s=size(state.imageProc.mathGUI.modImage);
	temp=reshape(state.imageProc.mathGUI.modImage,prod(s),1);
	[rows]=find(state.imageProc.mathGUI.modImage<0);
	temp(rows)=0;
	state.imageProc.mathGUI.modImage=reshape(temp,s(1),s(2));
end
loadImageFromArray('state.imageProc.mathGUI.modImage');
set(gh.mathGUI.figure1,'Pointer','arrow');
