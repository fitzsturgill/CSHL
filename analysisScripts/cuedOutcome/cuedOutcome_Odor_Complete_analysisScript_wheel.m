% cued outcome add wheel and pupil

all_behavior_trials = punishTrials;

    savename = 'allBehavior_rasters';
    ensureFigure(savename, 1);
%     reversals = find(diff(TE.BlockNumber(csPlusTrials, :))) + 1;
    
    subplot(1,3,1); 
    eventRasterFromTE(TE, all_behavior_trials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('Licks'); ylabel('trial number');
    set(gca, 'XLim', [-4 7]);
    
    
%     set(gca, 'YLim', [0 max(trialCount)]);
%     set(gca, 'FontSize', 14)
    wp = [bpX2pnt(-4, 20, -4) bpX2pnt(3, 20, -4)];
    subplot(1,3,2);    
    image(TE.Wheel.data.V(all_behavior_trials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
    set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 2, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 2]); 
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('Velocity');
    set(gca, 'YTick', []);
    
        subplot(1,3,3); phRasterFromTE(TE, all_behavior_trials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    set(gca, 'YTick', [], 'XLim', [-4 7]); 
    