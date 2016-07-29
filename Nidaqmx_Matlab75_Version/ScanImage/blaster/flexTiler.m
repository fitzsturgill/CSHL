function flexTiler(num_pixels, th_lo, th_hi)

global state;

state.blaster.flexTilerIm=zeros(num_pixels);

x_pixratio=size(state.acq.trackerReferenceAll, 1)/num_pixels;
y_pixratio=size(state.acq.trackerReferenceAll, 2)/num_pixels;

for i=1:num_pixels
    for j=1:num_pixels
        begin_x=round(1+(i-1)*x_pixratio);
        end_x=round(i*x_pixratio);
        begin_y=round(1+(j-1)*y_pixratio);
        end_y=round(j*y_pixratio);
        
        state.blaster.flexTilerIm(i,j)=mean2(state.acq.trackerReferenceAll(begin_x:end_x, begin_y:end_y));
    end
end

poscount=0;

for i=1:num_pixels
    for j=1:num_pixels
        if(state.blaster.flexTilerIm(i,j)>th_lo & state.blaster.flexTilerIm(i,j)<th_hi)
            % place blaster location here
            poscount=poscount+1;
            addNewBlasterPosByScreenXY(1+(j-1/2)*y_pixratio, 1+(i-1/2)*x_pixratio);
            %disp(['adding blaster location: ' num2str(i) ' ' num2str(j)])
            addEntryToNotebook(2, ['flextiler Pos ' num2str(poscount) ' (' ... 
            num2str(1+(i-1/2)*x_pixratio) ', ' num2str(1+(j-1/2)*y_pixratio) ')']);

        end
    end
end

setupBlasterConfigsByPos(22)
makeBlasterConfigMenu
updatereferenceimage