function setPath_dGK(incNum)
% function FSSetPath_dGK(incnum)
    global state
    
    if nargin < 1
        incNum=1;
    end
    
    basedir='N:\Microscope\Wenjia\dGK';
    cd(basedir);
    
    %today=datestr(now, 'ddmmmyy');
    
    notFoundDir=1;
    
    bool=1;
    
    while(notFoundDir)
       dirname=[basedir 'dGK_' num2str(incNum)];
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
               disp(['***Experiment Number too high. --- dGK_' num2str(incNum) ' does not exist on Microscope.--- Set to lower value.****']);
               return
           end
       end 
    end
    
    notFoundDir=1;
    basedir='C:\WY\dGK\';
    
    cd(basedir);
    
    while(notFoundDir)
       dirname=[basedir 'dGK_' num2str(incNum)];
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
    
    state.files.baseName=['dGK_' num2str(incNum) '_'];
    updateGuiByGlobal('state.files.baseName');
    
    disp(['*** SAVE PATH = ' state.files.savePath ' ***']);
    disp(['*** BASE NAME = ' state.files.baseName ' ***']);






