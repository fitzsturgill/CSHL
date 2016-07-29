function FSSetPath_AMPK(incNum)
% function FSSetPath_AMPK(incnum)
    global state
    
    if nargin < 1
        incNum=1;
    end
    
    basedir='N:\Microscope\Fitz\AMPK\';
    cd(basedir);
    
    %today=datestr(now, 'ddmmmyy');
    
    notFoundDir=1;
    
    bool=1;
    
    while(notFoundDir)
       dirname=[basedir 'AMPK_' num2str(incNum)];
	   if incNum ==1;
		   break
	   end
       try
           cd(dirname);
           incNum=incNum+1; %this will only happen if cd wass successful
           bool=0;
       catch
           notFoundDir=0;
           if bool
               disp(['***Experiment Number too high. --- AMPK_' num2str(incNum) ' does not exist on Microscope.--- Set to lower value.****']);
               return
           end
       end 
    end
    
    notFoundDir=1;
    basedir='C:\Fitz\AMPK\';
    
    cd(basedir);
    
    while(notFoundDir)
       dirname=[basedir 'AMPK_' num2str(incNum)];
       try
           cd(dirname);
           incNum=incNum+1; %this will only happen if cd wass successful
       catch
           notFoundDir=0;
       end 
    end    
    
    state.files.savePath=[dirname '\'];
    
    evalin('base', ['!mkdir ' state.files.savePath]);
    
    state.files.savePath=[dirname '\'];
	updateFullFileName(0);
	cd(state.files.savePath);
    
    state.files.baseName=['AMPK_' num2str(incNum) '_'];
    updateGuiByGlobal('state.files.baseName');
    
    disp(['*** SAVE PATH = ' state.files.savePath ' ***']);
    disp(['*** BASE NAME = ' state.files.baseName ' ***']);






