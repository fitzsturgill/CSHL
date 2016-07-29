function looppsd(num)

global state
    %kill('avPSD')
    
    waveo('avPSD', []);
    
    
    setWaveUserDataField('avPSD', 'nComponents', 0);
    %waveo('avPSD', []);
for i=1:num
   executegrabonecallback;
   waittilstop(state.daq.grabInput, 20);
   makepsd;
   avgin('PSD', 'avPSD');
end