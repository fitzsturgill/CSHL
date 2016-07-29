function plotMUAAverages(cycles, mode, channel, epoch)
    %cycles in rows are plotted together, those in columns, seperately
    % creates a new plot for the current epoch
    
    
    global state
    

    if nargin < 2
        mode = 'tile';
    end
    
    if nargin < 3
        channel = 'All';
    end
    
    if nargin < 4
        epoch = state.epoch;
    end
    
    lineColor = {[0 0 0] [1 0 0] [0 1 0] [0 1 1]  [1 0 1] [0 0 1]...
        [0 0 .5] [0 .5 0] [.5 0 0] [1 1 0] [0 0 0] [1 0 0] [0 1 0] [0 1 1]  [1 0 1] [0 0 1]};
    
    fig = figure('Name', ['MUA Averages: Epoch' num2str(epoch)]);
    
    for i = 1:size(cycles, 1)
        ax = axes('Parent', fig);
        title(['Cyc' num2str(cycles(i, :))]);
%         annotation('textbox', [0.1 0.7 0.5 0.2], 'string', ['Cyc' num2str(cycles(i, :))]);
        for j = 1:size(cycles, 2)
            cycle = cycles(i, j);
            if cycle == 0
                continue
            end
            wn = MUAAvgName(channel, cycle, epoch);
            if ~iswave(wn)
                newAverage(wn)
            end
            
            append(wn, 'axes', ax);
            setPlotProps(wn, 'Color', lineColor{1, j});
        end
    end
    
    switch lower(mode)
        case 'vertical'
            splayAxisVertical;
        case 'horizontal'
            splayAxisHorizontal;
        case 'tile'
            splayAxisTile;
        otherwise
            disp('Error in plotMUAAverages: incorect mode, mode must be "vertical", "horizontal", or "tile"');
            disp('Proceeding with Vertical Axis Splay Mode');
    end
            
            
            
            
            
            