%{
Any additional TE fields to add to cued outcome experiments should be added
in this script to keep everything organized so I don't have to sort through
all of my graphics code to find where I did actual analysis
%}

% settings for deconvolution
tau = 1; % time constant of exponential decay kernel for 6f see https://www.nature.com/articles/nature12354.pdf
kDuration = 3; % duration of kernel
Fs = 20;
kt = (1:kDuration*20)*(1/Fs) - (1/Fs); 
k = exp(-1 * (1/tau) * kt);
k = k / trapz(k);
epsilon = 0.1;
deconvSettings.epsilon = epsilon;
deconvSettings.tau = tau;
deconvSettings.kernelLength = kDuration;

DB = dbLoadExperiment('cuedOutcome');

for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    success = dbLoadAnimal(DB, animal);
    display(animal);
    nTrials = length(TE.filename);
    if ~success
        disp('wtf');
        continue
    end
    
    csWindow = cellfun(@(x,y) y(end) - x(1), TE.Cue, TE.Delay);
    csWindow = [zeros(nTrials, 1) csWindow];
    
    TE.lickLatency_cs = calcEventLatency(TE, 'Port1In', TE.Cue, TE.Us);    
    % deconvolution    
    TE.Photometry.data(1).dFdeconv = bpDeconv(TE.Photometry.data(1).raw - TE.Photometry.data(1).blF, k, epsilon, 'none');
    TE.Photometry.data(1).ZSdeconv = TE.Photometry.data(1).dFdeconv ./ nanmean(nanstd(TE.Photometry.data(1).dFdeconv(:,bpX2pnt(1, 20):bpX2pnt(4, 20)), 0, 2));
    TE.Photometry.settings.deconvSettings = deconvSettings;
    dbSaveAnimal(DB, animal); 
    

end