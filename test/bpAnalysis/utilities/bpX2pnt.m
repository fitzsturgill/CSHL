function p=bpX2pnt(x, Fs, startX)

    if nargin < 3
        startX = 0;
    end
    
    deltaX=1/Fs;


    p=max(round(1+(x-startX)/deltaX), 1);