function makepsd

global c2r1
global PSD

data=c2r1.data;

np=size(c2r1.data,2);

data=data(20:size(c2r1.data,2));

avg=mean(data);
data=data-avg;


Y=fft(data',np);
PSDm=Y.*conj(Y)/np;
PSDm=PSDm(1:(np/2+1));

PSDm=smooth2(PSDm,10);
waveo('PSD', PSDm);
PSD.xscale=[0 500/np];

