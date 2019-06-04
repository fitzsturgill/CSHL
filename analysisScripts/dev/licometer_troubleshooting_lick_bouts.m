    saveOn = 1;
    DB = dbLoadExperiment('reversals_noPunish_publish');

    nAnimals = length(DB.animals);
    nc = ceil(sqrt(nAnimals));
    nr = ceil(nAnimals/nc);

    savename = 'lick_bouts_reward_histograms';
    ensureFigure(savename, 1);
    ia2 = [];
    for counter = 1:length(DB.animals)
        animal = DB.animals{counter}
        dbLoadAnimal(DB, animal);
        goodTrials = find(rewardTrials);
        intervals = eventIntervalsFromTE(TE, 'Port1In', [0 1], TE.Us);

        ia = [];
        for counter2 = 1:length(goodTrials) 
            trial = goodTrials(counter2);
            ia = [ia intervals.intervals{trial}]; 
        end
        ia2 = [ia2 ia]; 
        subplot(nc,nr,counter); hold on;
        h = histogram(ia, 0:0.005:0.3);
    %     f = fit(h.BinEdges(1:end-1)' + h.BinWidth/2, h.BinCounts', 'gauss2');
    %     plot(f);
        textBox(animal);

        drawnow;
    end

    subplot(nc, nr, counter + 1); hold on;
    h=histogram(ia2, 0:0.005:0.3);

    f = fit(h.BinEdges(1:end-1)' + h.BinWidth/2, h.BinCounts', 'gauss2');
    rewLickRate = f.b2;
    plot(f);
    grid on;
    textBox(sprintf('%s, %.3g', 'all', f.b2));
    line([f.b2 f.b2], get(gca, 'YLim'));

