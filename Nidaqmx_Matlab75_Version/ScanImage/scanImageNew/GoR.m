function gr=GoR

global state
m=mean2(state.acq.acquiredData{1});
m0=0;
eval(['m0=state.acq.binFactor*state.acq.pmtOffsetChannel' num2str(1) ';']);
g=m-m0;

m=mean2(state.acq.acquiredData{2});
m0=0;
eval(['m0=state.acq.binFactor*state.acq.pmtOffsetChannel' num2str(2) ';']);
r=m-m0;

gr=g/r;

