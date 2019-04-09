%{
Any additional TE fields to add to cued outcome experiments should be added
in this script to keep everything organized so I don't have to sort through
all of my graphics code to find where I did actual analysis
%}



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
    
    dbSaveAnimal(DB, animal);      
end