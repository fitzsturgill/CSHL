function histAndNoiseAnalysisPixel
% This function uses the current ROI to determine the pixel by pixel mean and threshold image
% This will use the SEM which is the sd/sqrt(n)

global state gh
value = get(gh.imageProcessingGUI.fileName, 'Value');
[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
img=state.imageProc.currentImage(y:y1,x:x1,state.imageProc.currentFrame:state.imageProc.numberOfFrames);
n=1+state.imageProc.numberOfFrames-state.imageProc.currentFrame; 
if ~strcmp(class(img),'double')    % Not a double...
    class=class(img);
    img=double(img);
else
    class='double';
end

meanImg=mean(img,3);
stdImg=std(img,0,3);
stdImg=stdImg/sqrt(n);
thresholdImg= meanImg + state.imageProc.parsing.threshold*stdImg;	%This is the pixel by pixel threshold
sizeImg=size(thresholdImg);
pixels=sizeImg(1)*sizeImg(2);
rowOfImg=reshape(thresholdImg,pixels,1);   % MAde a columns vector
rowOfImg(isnan(rowOfImg))=0;	%Remove NANs
eval(['rowOfImg = ' class '(rowOfImg);']);
state.imageProc.parsing.thresholdImg=rowOfImg;
rowOfImg = reshape(rowOfImg,sizeImg(1),sizeImg(2));

%sizeImg=size(meanImg);
%pixels=sizeImg(1)*sizeImg(2);
%meanImg=double(meanImg);
rowOfmeanImg=reshape(meanImg,pixels,1);   % MAde a columns vector
rowOfmeanImg(isnan(rowOfmeanImg))=0;	%Remove NANs
eval(['rowOfmeanImg = ' class '(rowOfmeanImg);']);
state.meanImg=rowOfmeanImg;

%figure('DoubleBuffer','On','NumberTitle','off','Name','Pixel Histogram and Image with Noise Removed',...
%    'Color','white','pos',[106   507   868   389]);
%cmap = [ .9 0 .9;gray(255)]; 
%colormap(cmap);
%imagesc(medfilt2(rowOfImg,[3 3]));
%imagesc(rowOfImg);
%title(['Filtered Baseline Image with Intensities > ' num2str(state.imageProc.parsing.threshold) '  SEM above Mean'],'FontName','TimesNewRoman');
disp(['Created a baseline image with a threshold of mean + ' num2str(state.imageProc.parsing.threshold) '*SEM.']);
