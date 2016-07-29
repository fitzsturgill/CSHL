
function switcher(Tags)

global head_info target_directory extension filename chan imbuffer  
global anabuffer projRange
global contrastRange

switch Tags(1:5)   
   case 'Sweep' %if need to get new file:
   	SweepHandle=findobj(gcbf, 'Tag', 'SweepText'); 
      sweepStr=get(SweepHandle, 'String');
      sweep=str2num(sweepStr);   %note: need the transpose here
      DirectoryHandle=findobj(gcbf, 'Tag', 'DirectoryText');     
		directory=get(DirectoryHandle, 'String');                    
		BaseHandle=findobj(gcbf, 'Tag', 'BaseText');    
      base=get(BaseHandle, 'String');
      ChanHandle=findobj(gcbf, 'Tag', 'ChanText'); 
      chanStr=get(ChanHandle, 'String');
      chan=str2num(chanStr);
      if Tags == 'SweepNext'
         sweep=GetNextSweep(directory, base, sweep, extension, 1);
      	sweepStr=num2str(sweep);
         set(SweepHandle, 'String' ,sweepStr);
      end
      if Tags == 'SweepLast'
         sweep=GetNextSweep(directory, base, sweep, extension, -1);
      	sweepStr=num2str(sweep);
         set(SweepHandle, 'String' ,sweepStr);
      end
      %we migth be able to move this out of the case statement
   	filename=sprintf('%s%03d.cfd',strcat(directory,base),sweep);
   	head_info=cf4header(filename);
      imbuffer=CFDread2(filename,1,chan);  
      %max(squeeze(reshape(imbuffer,512*512,1)))
      %whos imbuffer
      imshow(imbuffer, contrastRange); 

      %truesize([512,512]);
      %setup the slider
		FrameSliderHandle=findobj(gcbf, 'Tag', 'FrameSlider');
		set(FrameSliderHandle, 'SliderStep',[1/(head_info.n_images-1),10/(head_info.n_images-1)]); 
		FrameMaxHandle=findobj(gcbf, 'Tag', 'FrameMax');
		set(FrameMaxHandle, 'String', num2str(head_info.n_images));
		FrameMinHandle=findobj(gcbf, 'Tag', 'FrameMin');
		set(FrameMinHandle, 'String', '1');
		FrameCurrentHandle=findobj(gcbf, 'Tag', 'FrameCurrent');
		set(FrameCurrentHandle, 'String', '1');
      
   case 'Frame'%if need to get new frame in the same file:
        if Tags(1:11) == 'FrameCurren'
         FrameCurrentHandle=findobj(gcbf, 'Tag', 'FrameCurrent');
         frameStr=get(FrameCurrentHandle, 'String');
         Frame=str2num(frameStr);
         FrameSliderHandle=findobj(gcbf, 'Tag', 'FrameSlider');
	      set(FrameSliderHandle, 'Value',(Frame-1)/(head_info.n_images-1));
      end
      if Tags(1:11) == 'FrameSlider';
         FrameSliderHandle=findobj(gcbf, 'Tag', 'FrameSlider');
         value=get(FrameSliderHandle, 'Value');
         Frame=round(value*(head_info.n_images-1)+1);
         FrameCurrentHandle=findobj(gcbf, 'Tag', 'FrameCurrent');
         set(FrameCurrentHandle, 'String', num2str(Frame));
      end
      if Tags(1:11) == 'FrameMovies';
         ToggleHandle=findobj(gcbf, 'Tag', 'FrameMovies');
         toggleValue=get(ToggleHandle, 'value');
         MovieTimingHandle=findobj(gcbf, 'Tag', 'FrameFPS');
         valueStr=get(MovieTimingHandle, 'String');
			value=str2num(valueStr);
			interval=1/value;
         FrameSliderHandle=findobj(gcbf, 'Tag', 'FrameSlider');
         FrameCurrentHandle=findobj(gcbf, 'Tag', 'FrameCurrent');
         i=1;
         imbuffer=CFDread2(filename,i,chan);  
         h=imshow(imbuffer, contrastRange);
         set(h,'erasemode','xor'); % to avoid image flicker
         i=i+1;
         while (i < (head_info.n_images+1)) & (toggleValue == 1)
         	value=(i-1)/(head_info.n_images-1);
            set(FrameSliderHandle, 'Value', value);
            value=num2str(i);
            set(FrameCurrentHandle, 'String', value);
            imbuffer=CFDread2(filename,i,chan);  
            set(h, 'cdata', imbuffer);  % use this instead of imshow for flickerless display
            %imshow(imbuffer, []);
            pause(interval);
            i=i+1;
            toggleValue=get(ToggleHandle, 'value'); % to turn off movie with button
         end
         Frame=1;
      end
      
      if Tags(1:11) == 'FrameFocusT';
         for i=1:10; 
	         imbuffer=CFDread2(filename,i*4,chan);  
            imshow(imbuffer, contrastRange); 
            pause(0);
            Frame=1;
			end				         
      end
      
      imbuffer=CFDread2(filename,Frame,chan);  
      imshow(imbuffer, contrastRange);                           
      %     truesize( [512,512]);   
   case 'Iprop' %display and image properties
      switch Tags
      case 'IpropMaxSl'
         sliderHandle=findobj(gcbf, 'Tag', Tags);
         value=get(sliderHandle, 'Value');
         contrastRange(2)=round((value*255)+1);
         inHandle=findobj(gcbf, 'Tag','IpropMaxIn');
         valueStr=num2str(contrastRange(2));
         set(inHandle, 'String', valueStr);
         % 'here'
          contrastRange
      case 'IpropMinSl'
         sliderHandle=findobj(gcbf,'Tag', Tags);
         value=get(sliderHandle, 'Value');
         contrastRange(1)=round((value*255)+1);
         if contrastRange(1) > (contrastRange(2)-1) 
            contrastRange(1) = contrastRange(2)-1;
            set(sliderHandle, 'Value', (contrastRange(1)-1)/255);
         end
         inHandle=findobj(gcbf, 'Tag','IpropMinIn');
         valueStr=num2str(contrastRange(1));
         set(inHandle, 'String', valueStr);
         'here'
          contrastRange
      otherwise
      end
      imshow(imbuffer, contrastRange);
               
				      
   case 'Proje' %if need to make projections
      ChanHandle=findobj(gcbf, 'Tag', 'ChanText'); 
      chanStr=get(ChanHandle, 'String');
      chan=str2num(chanStr);
      pointerStr=get(gcf, 'Pointer');
      set(gcf, 'Pointer', 'watch');
      switch Tags
      case 'ProjeZ'
         imbuffer=ProjStack(filename, chan, projRange, 'z');
         whos imbuffer;
      case 'ProjeX'
         imbuffer=ProjStack(filename, chan, projRange, 'x');
      case 'ProjeY'
         imbuffer=ProjStack(filename, chan, projRange, 'y');
      otherwise 
      end
      imshow(imbuffer, contrastRange);  
      %'yes'
      %      truesize( [512,512]);
      set(gcf, 'Pointer', 'default');
		      
   case 'Saves' % if need to save data currently on the screen
      switch Tags
      case 'SavesImg'
         outname=NextFilename(target_directory, filename, '.tif')
         imwrite(imbuffer,outname,'tif');
      case 'SavesDat'
         outname=NextFilename(target_directory, filename, '.ana');
         fid = fopen(outname,'w');
         fprintf(fid,'%6.2f\n',anabuffer);
         fclose(fid);
      otherwise
      end
      
   case 'Analy' % for analysis
      %DatHandle=findobj('Tag', 'WinDat');
      %axes(DatHandle);
      switch Tags % put analysis code here
      case 'AnalyBox'
         anabuffer=boxScan(filename);
		case 'AnalyLin'
         anabuffer=lineScan(filename);
      otherwise
      end
         DatHandle=findobj('Tag', 'WinDat');
      	axes(DatHandle);
         plot(anabuffer, 'Color', 'y');
         ImgHandle=findobj('Tag', 'WinImg');
         axes(ImgHandle);
         % new stuff
      case 'Anima' % for analysis
         switch Tags % put analysis code here
         case 'AnimaAvi'
            pointerStr=get(gcf, 'Pointer');
		      set(gcf, 'Pointer', 'watch');
            MovieTimingHandle=findobj(gcbf, 'Tag', 'FrameFPS');
         	valueStr=get(MovieTimingHandle, 'String');
				rate=str2num(valueStr);
            WriteMovie(filename, rate, '.avi');
            set(gcf, 'Pointer', pointerStr);
         otherwise
         end
         
         %anabuffer=boxScan(filename);

         
      otherwise      
      end
      
% ------------------------------------------------------------------------------
% local functions
% ------------------------------------------------------------------------------

function yy=WriteMovie(filename, rate, type)
	
      global head_info target_directory extension filename chan imbuffer  
      
      xx(1:head_info.pixels_xy(1), 1:head_info.pixels_xy(2), 1:head_info.n_images)=uint8(0);
      xx=uint8(xx);
      for i=1:head_info.n_images; xx(:, :, i)= CFDread2(filename, i, chan); end
      fn_out=NextFilename(target_directory, filename, type);
      AviWrite(xx,fn_out,rate);
      yy=1;
    