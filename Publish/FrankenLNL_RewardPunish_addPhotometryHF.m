% add PhotometryHF field to TE that isn't downsampled so severely

DB = dbLoadExperiment('FrankenLNL_RewardPunish');

BL = [1 4];
dFFMode = 'expFit'; % for both channels
for counter = length(DB.animals) - 1:length(DB.animals)
    animal = DB.animals{counter}
    dbLoadAnimal(DB, animal);    
    sessions = bpLoadSessionsFromTE(TE);
    channels = sessions(1).SessionData.NidaqData{1,2}.channelsOn;   
    TE.PhotometryHF = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue2', 'channels', channels, 'baseline', BL, 'downsample', 61);
    dbSaveAnimal(DB, animal);            
end