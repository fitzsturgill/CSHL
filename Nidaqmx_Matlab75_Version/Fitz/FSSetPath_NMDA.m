function FSSetPath_phy(incNum)
% function FSSetPath_phy(incnum)
    global state
    
    if nargin < 1
        incNum=1;
    end
    
    basedir='N:\Microscope\Fitz\phy\';
    cd(basedir);
    
    %today=datestr(now, 'ddmmmyy');
    
    notFoundDir=1;
    
    bool=1;
    
    while(notFoundDir)
       dirname=[basedir 'phy_' num2str(incNum)];
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
               disp(['***Experiment Number too high. --- phy_' num2str(incNum) ' does not exist on Microscope.--- Set to lower value.****']);
               return
           end
       end 
    end
    
    notFoundDir=1;
    basedir='C:\Fitz\phy\';
    
    cd(basedir);
    
    while(notFoundDir)
       dirname=[basedir 'phy_' num2str(incNum)];
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
    
    state.files.baseName=['phy_' num2str(incNum) '_'];
    updateGuiByGlobal('state.files.baseName');
    
    disp(['*** SAVE PATH = ' state.files.savePath ' ***']);
    disp(['*** BASE NAME = ' state.files.baseName ' ***']);






% function FSSetPath_mGluR
% 
%     global state
%     
%     basedir='C:\Fitz\mGluR\';
%     cd(basedir);
%     
%     %today=datestr(now, 'ddmmmyy');
%     
%     notFoundDir=1;
%     incNum=76;
%     
%     while(notFoundDir)
%        dirname=[basedir 'mGluR_' num2str(incNum)];
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
%     state.files.baseName=['mGluR_' num2str(incNum) '_'];
%     updateGuiByGlobal('state.files.baseName');
%     
%     disp(['*** SAVE PATH = ' state.files.savePath ' ***']);
%     disp(['*** BASE NAME = ' state.files.baseName ' ***']);