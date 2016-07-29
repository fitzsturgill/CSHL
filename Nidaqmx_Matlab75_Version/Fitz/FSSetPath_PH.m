function FSSetPath(incNum)
% function FSSetPath(incnum)
    global state
    
    if nargin < 1
        incNum=100;
    end
    
    basedir='N:\Microscope\Fitz\PH\PH\';
    cd(basedir);
    
    %today=datestr(now, 'ddmmmyy');
    
    notFoundDir=1;
    
    bool=1;
    
    while(notFoundDir)
       dirname=[basedir 'PH_' num2str(incNum)];
       try
           cd(dirname);
           incNum=incNum+1; %this will only happen if cd wass successful
           bool=0;
       catch
           notFoundDir=0;
           if bool
               disp(['***Experiment Number too high. --- PH_' num2str(incNum) ' does not exist.--- Set to lower value.****']);
               return
           end
       end 
    end
    
    notFoundDir=1;
    basedir='C:\Fitz\PH\';
    
    cd(basedir);
    
    while(notFoundDir)
       dirname=[basedir 'PH_' num2str(incNum)];
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
    
    state.files.baseName=['PH_' num2str(incNum) '_'];
    updateGuiByGlobal('state.files.baseName');
    
    disp(['*** SAVE PATH = ' state.files.savePath ' ***']);
    disp(['*** BASE NAME = ' state.files.baseName ' ***']);


% function FSSetPath
% 
%     global state
%     
%     basedir='C:\Fitz\PH\';
%     cd(basedir);
%     
%     %today=datestr(now, 'ddmmmyy');
%     
%     notFoundDir=1;
%     incNum=219;
%     
%     while(notFoundDir)
%        dirname=[basedir 'PH_' num2str(incNum)];
%        try
%            cd(dirname);
%            incNum=incNum+1; %this will only happen if cd wass successful
%        catch
%            notFoundDir=0;
%        end 
%     end
%     
%     state.files.savePath=[dirname '\'];
%     
%     evalin('base', ['!mkdir ' state.files.savePath]);
%     
%     state.files.savePath=[dirname '\'];
% 	updateFullFileName(0);
% 	cd(state.files.savePath);
%     
%     state.files.baseName=['PH_' num2str(incNum) '_'];
%     updateGuiByGlobal('state.files.baseName');
%     
%     disp(['*** SAVE PATH = ' state.files.savePath ' ***']);
%     disp(['*** BASE NAME = ' state.files.baseName ' ***']);