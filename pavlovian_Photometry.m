function pavlovian_Photometry(SessionData, sessionTitle)
    fig = figure('Name', sessionTitle);
    title_ax = textAxes(fig, sessionTitle, 12); %fontsize = 12
    
    stimBarY = -0.05;
    yLims = [-0.1 0.1];

    %primary comparison: toneA reward vs toneB punish

    ax(1) = axes; hold on
    sessionAvg(SessionData, 1, 2, 'g', ax(1));
    sessionAvg(SessionData, 2, 4, 'r', ax(1));
    textBox('green: toneA, reward; red: toneB, punish', [],[], 12);
    xlabel('time (s)');
    ylabel('GCaMP norm');
    addStimulusBar(ax(1), [2 3 stimBarY], 'tone', [0 0 0], 1);
    addStimulusBar(ax(1), [3.5 3.7 stimBarY], 'feedback', [0 0 0], 1);


    % toneA: reward, ommision, punish
    ax(2) = axes; hold on  

    sessionAvg(SessionData, 1, 3, 'k', ax(2));
    sessionAvg(SessionData, 1, 4, 'r', ax(2));
    sessionAvg(SessionData, 1, 2, 'g', ax(2));
    textBox('Tone A', [],[], 12);
    xlabel('time (s)');
    ylabel('GCaMP norm');
    addStimulusBar(ax(2), [2 3 stimBarY], 'tone', [0 0 0], 3);
    addStimulusBar(ax(2), [3.5 3.7 stimBarY], 'feedback', [0 0 0], 3);

    % toneB: reward, ommision, punish
    ax(3) = axes; hold on 
    sessionAvg(SessionData, 2, 2, 'g', ax(3));
    sessionAvg(SessionData, 2, 3, 'k', ax(3));
    sessionAvg(SessionData, 2, 4, 'r', ax(3));
    textBox('Tone B', [],[], 12);
    xlabel('time (s)');
    ylabel('GCaMP norm');
    addStimulusBar(ax(3), [2 3 stimBarY], 'tone', [0 0 0], 3);
    addStimulusBar(ax(3), [3.5 3.7 stimBarY], 'feedback', [0 0 0], 3);
    
    set(ax, 'YLim', yLims);
    splayAxisTile;
    
    
    
    figure('Name', 'SupriseFigure');
    ax2(1) = axes; hold on;

    sessionAvg(SessionData, 1, 2, 'g', ax2(1));
    sessionAvg(SessionData, 2, 2, 'k', ax2(1));
    textBox('reward suprise (A vs B)', [],[], 12);
    xlabel('time (s)');
    ylabel('GCaMP norm');
    addStimulusBar(ax2(1), [2 3 stimBarY], 'tone', [0 0 0], 3);
    addStimulusBar(ax2(1), [3.5 3.7 stimBarY], 'feedback', [0 0 0], 3);
    
    
    ax2(2) = axes; hold on;

    sessionAvg(SessionData, 2, 4, 'r', ax2(2));
    sessionAvg(SessionData, 1, 4, 'k', ax2(2));
    textBox('punishment suprise (B vs A)', [],[], 12);
    xlabel('time (s)');
    ylabel('GCaMP norm');
    addStimulusBar(ax2(2), [2 3 stimBarY], 'tone', [0 0 0], 3);
    addStimulusBar(ax2(2), [3.5 3.7 stimBarY], 'feedback', [0 0 0], 3);
    set(ax2, 'YLim', yLims);
    splayAxisTile;
    
    
    
% fig = figure;
% lickRaster(SessionData, 1, 2, fig, 'toneA reward')
% lickRaster(SessionData, 1, 3, fig, 'toneA omission')
% lickRaster(SessionData, 1, 4, fig, 'toneA punish')
% 
% lickRaster(SessionData, 2, 2, fig, 'toneB reward')
% lickRaster(SessionData, 2, 3, fig, 'toneB omission')
% lickRaster(SessionData, 2, 4, fig, 'toneB punish')
% 
% splayAxisTile;