

trials = 1:length(TE.filename);
fh = zeros(length(trials), 1);
for counter = 1:length(trials)
    trial = trials(1);
    saveName = sprintf('wheel_example_traces_tr%d', trial);
    fh(counter) = ensureFigure(saveName, 1);
    mcLandscapeFigSetup(gcf);
    ax = axes; hold on;
    Fs = 20;
    duration = 120;

    data = zeros(Fs*duration, 5); % wheel pupil whisk chat dat
    data(:,1) = nanzscore(TE.Wheel.data.V(trial,1:Fs*duration));
    data(:,2) = nanzscore(TE.Whisk.whiskNorm(trial,1:Fs*duration));
    data(:,3) = nanzscore(TE.pupil.pupDiameterNorm(trial,1:Fs*duration));
    data(:,4) = nanzscore(TE.Photometry.data(1).ZS(trial,:));
    data(:,5) = nanzscore(TE.Photometry.data(2).ZS(trial,:));
    xdata = TE.Photometry.xData;

    % splay the data vertically
    splay = range(data);
    splay = cumsum(splay);
    data = data - splay;

    rewX = TE.Reward{trial}(:,1) - TE.Photometry.startTime(trial);
    rewX = repmat(rewX', 2, 1);
    rewY = repmat([max(data(:)); min(data(:))], 1, size(rewX, 2));

    plot(rewX, rewY, 'Color', [0.7 0.7 0.7]);
    lh = plot(xdata, data, 'k');
    lh(4).Color = mycolors('chat');
    lh(5).Color = mycolors('dat');
    ax.YAxis.Visible = 'off';
    ax.XLim = [0 120];
    xlabel('time (s)');
end
