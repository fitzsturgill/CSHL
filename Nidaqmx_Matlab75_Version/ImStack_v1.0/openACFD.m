function openACFD(filename)
global state gh head_info

head_info = CF4header(filename); 
state.imageProc.cfd.image = [];
totalimages=head_info.n_images;

state.imageProc.cfd.image = openACFDSpine(filename, state.imageProc.cfd.numberofChannels);
loadImageFromArray('state.imageProc.cfd.image');

