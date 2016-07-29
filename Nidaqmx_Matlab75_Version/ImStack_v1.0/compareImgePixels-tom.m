function compareImgePixels
% This is the function that will read in an image and dispaly a histogram of
% the pixels as well as the rescaled the histogram.
% Ti then looks for only thjose pixel intensities greater than 5 std above the mean and
% displays this image.
global state gh

% Load the dat from the software...
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

% Now we have the dta aavergaed over a few frames.  This is the signal.


% Display Mean Image
%figure('DoubleBuffer','On','NumberTitle','off','Name','Pixel Histogram and Image with Noise Removed',...
%    'Color','white','pos',[106   507   868   389]);
%cmap = [ .9 0 .9;gray(255)]; 
%colormap(cmap);
%eval(['meanImg = ' class '(meanImg);']);
% imagesc(medfilt2(rowOfImg,[3 3]));
%imagesc(medfilt2(meanImg,[3 3]));
%title(['Filtered Mean Image'],'FontName','TimesNewRoman');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rowOfImg=reshape(meanImg,pixels,1);   % MAde a columns vector
rowOfImg(isnan(rowOfImg))=0;	%Remove NANs
if length(rowOfImg) ~= length(state.imageProc.parsing.thresholdImg)
	disp('Must choose baseline image first');
	return
end

% Copy the data
copyrowOfImg=rowOfImg; 	%RowOfImge is the good data.

%Find out which elements of hte data match the criteria.
if state.imageProc.parsing.threshold > 0 
	r=find(rowOfImg < double(state.imageProc.parsing.thresholdImg));
else
	r=find(rowOfImg > double(state.imageProc.parsing.thresholdImg));
end

% The other ones....
if state.imageProc.parsing.threshold > 0 
	r2=find(rowOfImg >= double(state.imageProc.parsing.thresholdImg));
else
	r2=find(rowOfImg <= double(state.imageProc.parsing.thresholdImg));
end

%Set thos pixels to 0
rowOfImg(r) = 0;

% correct for ther mean value....
meanImg=double(state.meanImg);
meanImg(r2)=0;

rowOfImg=rowOfImg-double(state.meanImg);
r=find(rowOfImg<0);
rowOfImg(r) = 0;

r=find(copyrowOfImg >= double(state.imageProc.parsing.thresholdImg));
copyrowOfImg(r) = 0;
crowOfImg=rowOfImg;

% rowOfImg = reshape(rowOfImg,sizeImg(1),sizeImg(2));

% state.imageProc.internal.thresholdImage=rowOfImg;
% loadImageFromArray('state.imageProc.internal.thresholdImage');

% copyrowOfImg = reshape(copyrowOfImg,sizeImg(1),sizeImg(2));

figure('DoubleBuffer','On','NumberTitle','off','Name','Overlay: Signal over Mean Image.',...
   'Color','white','pos',[106   507   868   389]);
red=makeColorMap('red',7);
green=makeColorMap('green',7);
cmap = [red;green]; 
colormap(cmap);

minMean=min(meanImg);
meanImg=(meanImg-minMean);
maxMean=max(meanImg);
meanImg=(meanImg)*127/maxMean;

minSignal=min(crowOfImg);
meanSignal=(crowOfImg-minSignal);
maxSignal=max(meanSignal);
signal=(meanSignal)*127/maxMean + 128;

new=signal+meanImg;
new=reshape(new,sizeImg(1),sizeImg(2));
eval(['new = ' class '(new);']);
imagesc(new);


% figure('DoubleBuffer','On','NumberTitle','off','Name','Pixel Histogram and Image with Noise Removed',...
%    'Color','white','pos',[106   507   868   389]);
% cmap = [ .9 0 .9;gray(255)]; 
% colormap(cmap);
% imagesc(medfilt2(rowOfImg,[3 3]));
% if state.imageProc.parsing.threshold > 0 
% 	title(['Median Filtered Image with Intensities > ' num2str(abs(state.imageProc.parsing.threshold)) '  SEM above Mean'],'FontName','TimesNewRoman');
% else
% 	title(['Median Filtered Image with Intensities < ' num2str(abs(state.imageProc.parsing.threshold)) '  SEM below Mean'],'FontName','TimesNewRoman');
% end
% 
% 
% figure('DoubleBuffer','On','NumberTitle','off','Name','Pixel Histogram and Image with Noise Removed',...
%    'Color','white','pos',[106   507   868   389]);
% cmap = [ .9 0 .9;gray(255)]; 
% colormap(cmap);
% imagesc(medfilt2(copyrowOfImg,[3 3]));
% if state.imageProc.parsing.threshold > 0 
% 	title(['Median Filtered Image with Intensities < ' num2str(abs(state.imageProc.parsing.threshold)) '  SEM above Mean'],'FontName','TimesNewRoman');
% else
% 	title(['Median Filtered Image with Intensities > ' num2str(abs(state.imageProc.parsing.threshold)) '  SEM below Mean'],'FontName','TimesNewRoman');
% end








