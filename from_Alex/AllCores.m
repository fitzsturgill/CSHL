disp('Recruiting six additional workers to serve you, master. One moment, please.')
warning off
if matlabpool('size') == 0
    matlabpool local 7
end
warning on
disp('Master, *WE* are here to serve you.')