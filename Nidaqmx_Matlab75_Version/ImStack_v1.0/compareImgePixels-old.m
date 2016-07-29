function compareImgePixels
% This is the function that will read in an image and dispaly a histogram of
% the pixels as well as the rescaled the histogram.
% Ti then looks for only thjose pixel intensities greater than 5 std above the mean and
% displays this image.
global state gh

value = get(gh.imageProcessingGUI.fileName, 'Value');
[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
img=state.imageProc.currentImage(y:y1,x:x1,state.imageProc.currentFrame:state.imageProc.numberOfFrames);
if ~strcmp(class(img),'double')    % Not a double...
    class=class(img);
    img=double(img);
else
    class='double';
end
meanImg=mean(img,3);
sizeImg=size(meanImg);

if length(sizeImg) > 2
    error('processImageDataYalin: Only works on 2 D data');
end
pixels=sizeImg(1)*sizeImg(2);

% Display Mean Image
figure('DoubleBuffer','On','NumberTitle','off','Name','Pixel Histogram and Image with Noise Removed',...
    'Color','white','pos',[106   507   868   389]);
cmap = [ .9 0 .9;gray(255)]; 
colormap(cmap);
eval(['meanImg = ' class '(meanImg);']);
% imagesc(medfilt2(rowOfImg,[3 3]));
imagesc(medfilt2(meanImg,[3 3]));
title(['Filtered Mean Image'],'FontName','TimesNewRoman');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meanImg=double(meanImg);
rowOfImg=reshape(meanImg,pixels,1);   % MAde a columns vector
rowOfImg(isnan(rowOfImg))=0;	%Remove NANs

if length(rowOfImg) ~= length(state.imageProc.parsing.thresholdImg)
	disp('Must choose baseline image first');
	return
end

copyrowOfImg=rowOfImg;
if state.imageProc.parsing.threshold > 0 
	r=find(rowOfImg < double(state.imageProc.parsing.thresholdImg));
	r2=find(rowOfImg > double(state.imageProc.parsing.thresholdImg));
else
	r=find(rowOfImg > double(state.imageProc.parsing.thresholdImg));
	r2=find(rowOfImg < double(state.imageProc.parsing.thresholdImg));
end
rowOfImg(r) = 0;
% r=find(copyrowOfImg >= double(state.imageProc.parsing.thresholdImg));
copyrowOfImg(r2) = 0;

newCopy=copyrowOfImg;
newRowCopy=rowOfImg;

rowOfImg = reshape(rowOfImg,sizeImg(1),sizeImg(2));
copyrowOfImg = reshape(copyrowOfImg,sizeImg(1),sizeImg(2));
eval(['rowOfImg = ' class '(rowOfImg);']);
eval(['copyrowOfImg = ' class '(copyrowOfImg);']);

figure('DoubleBuffer','On','NumberTitle','off','Name','Pixel Histogram and Image with Noise Removed',...
    'Color','white','pos',[106   507   868   389]);
cmap = [ .9 0 .9;zeros(255,3)]; 
colormap(cmap);
%c=zeros(256,3);
%c(1,:) = [.9 0 .9];
%colormap(c);
imagesc(medfilt2(rowOfImg,[3 3]));

if state.imageProc.parsing.threshold > 0 
	title(['Median Filtered Image with Intensities > ' num2str(abs(state.imageProc.parsing.threshold)) '  SEM above Mean'],'FontName','TimesNewRoman');
else
	title(['Median Filtered Image with Intensities < ' num2str(abs(state.imageProc.parsing.threshold)) '  SEM below Mean'],'FontName','TimesNewRoman');
end


% figure('DoubleBuffer','On','NumberTitle','off','Name','Pixel Histogram and Image with Noise Removed',...
%     'Color','white','pos',[106   507   868   389]);
% cmap = [ .9 0 .9;gray(255)]; 
% colormap(cmap);
% imagesc(medfilt2(copyrowOfImg,[3 3]));
% if state.imageProc.parsing.threshold > 0 
% 	title(['Median Filtered Image with Intensities < ' num2str(abs(state.imageProc.parsing.threshold)) '  SEM above Mean'],'FontName','TimesNewRoman');
% else
% 	title(['Median Filtered Image with Intensities > ' num2str(abs(state.imageProc.parsing.threshold)) '  SEM below Mean'],'FontName','TimesNewRoman');
% end

figure('DoubleBuffer','On','NumberTitle','off','Name','Pixel Histogram and Image with Noise Removed',...
    'Color','white','pos',[106   507   868   389]);
red=makeColorMap('red',7);
green=makeColorMap('green',7);
cmap=[red;green];
colormap(cmap);

newCopy=double(newCopy);
newRowCopy=double(newRowCopy);

maxP=max(max(newRowCopy));
minP=min(min(newRowCopy));
newRowCopy=(newRowCopy-minP)*128/maxP;

maxP=max(max(newCopy));
minP=min(min(newCopy));
newCopy=(newCopy-minP)*128/maxP + 129;

new=newCopy+newRowCopy;
% subplot(2,1,1);
new = reshape(new,sizeImg(1),sizeImg(2));
% new=medfilt2(new,[3 3]);
imagesc(new);
% subplot(2,1,2);
% imagesc(medfilt2(rowOfImg,[3 3]));
title('New Overlay for Yalin');









