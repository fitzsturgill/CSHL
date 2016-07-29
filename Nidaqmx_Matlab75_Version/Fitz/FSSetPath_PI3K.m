function FSSetPath

    global state
    
    basedir='C:\Fitz\PI3K\';
    cd(basedir);
    
    %today=datestr(now, 'ddmmmyy');
    
    notFoundDir=1;
    incNum=1;
    
    while(notFoundDir)
       dirname=[basedir 'PI3K_' num2str(incNum)];
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
    
    state.files.baseName=['PI3K_' num2str(incNum) '_'];
    updateGuiByGlobal('state.files.baseName');
    
    disp(['*** SAVE PATH = ' state.files.savePath ' ***']);
    disp(['*** BASE NAME = ' state.files.baseName ' ***']);