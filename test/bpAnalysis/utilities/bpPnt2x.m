function x=bpPnt2x(p, Fs, startX)
    
    if nargin < 3
        startX = 0;
    end
    deltaX=1/Fs;    
    x = startX + deltaX * (p - 1);
	