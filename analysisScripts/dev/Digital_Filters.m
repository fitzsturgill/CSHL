    %% try filtering first
    % note-   5 pole Butterworth filter in Matlab used in Frohlich and McCormick  
     % Create butterworth filter
    lowCutoff = 15/6100 * 2; % multiply by 2 to convert to rad/sample- see butter documentation
    % for a cutoff freq of 300Hz and sample rate of 1000Hz, cutoff
    % corresponds to 0.6pi rad/sample    300/1000 * 2 = 0.6    
    [b, a] = butter(10, lowCutoff, 'low'); % mimic filt filt in which filter order is doubled...
    
    [z,p,k] = butter(10,lowCutoff, 'low');
    sos = zp2sos(z,p,k);
    ensureFigure('test', 1);    
    hfvt = fvtool(b,a,sos);
    legend(hfvt, 'TF Design', 'ZPK Design');
    
    %%
    [b2, a2] = cheby1(10,3,lowCutoff);
    
    
    
    
%%
    
    
    