function phSessionRaster(analysis, type, outcome, fCh, varargin)

% [ax, hl] = phPlotSessionAverage(analysis, type, outcome, fCh, varargin);

    
    %% default values
    
    ax = [];
    fig = gcf;
    trials = -1; % -1 = no trial filtering, default
    secondsPerTick = 2; %
    downsample = 1; % factor by which to decimate/downsample
    lookupFactor = 4; % multiple of SD above and below mean to set CLim
    CLim = []; % exclusive with lookupFactor, manually supply CLim for image
    % parse input parameter pairs
    counter = 1;
    while counter+1 <= length(varargin) 
        prop = varargin{counter};
        val = varargin{counter+1};
        switch prop
            case 'linespec'
                linespec = val;
                if ischar(linespec)
                    linespec = {linespec};
                end
            case 'ax'
                ax = val;
            case 'fig'
                fig = val;
            case 'secondsPerTick'
                secondsPerTick = val;
            case 'downsample'
                downsample = val;
            case 'lookupFactor'
                lookupFactor = val;
            case 'CLim'
                CLim = val;
            otherwise
        end
        counter=counter+2;
    end

    %%
    
    if isempty(ax) %make a new axes unless one is provided
        ax=axes(...
        'Parent', fig,...
        'YDir', 'Reverse',...
        'TickDir', 'out'...
        );
    else
        set(ax,...
        'YDir', 'Reverse',...
        'TickDir', 'out'...
        );
        axes(ax);
    end
    

    if trials == -1
        cData = analysis.Photometry.data(type, outcome).data(:,:,fCh);
    else
        cData = analysis.Photometry.data(type, outcome).avg(trials,:,fCh);        
    end
    
    if downsample > 1
        for counter = 1:size(cData, 1)
            dd = decimate(cData(counter,:), downsample);
            if counter == 1
                cCopy = zeros(size(cData, 1), size(dd, 2));
            end
            cCopy(counter,:) = decimate(cData(counter,:), downsample);            
        end
        cData = cCopy;
    end
        
        
    x1 = min(analysis.Photometry.data(type, outcome).x(:,1));
    x2 = max(analysis.Photometry.data(type, outcome).x(:,1));
    y1 = 1;
    y2 = size(cData, 1);
    
    ih = image('XData', [x1 x2], 'YData', [y1 y2], 'CData', cData, 'CDataMapping', 'Scaled', 'Parent', ax);
    set(ax, 'XLim', [x1 x2], 'YLim', [y1 y2]);
    
    if ~isempty(lookupFactor)
        m = mean(mean(cData, 'omitnan'));
        s = mean(std(cData, 'omitnan'));
        set(ax, 'CLim', [m - lookupFactor * s, m + lookupFactor * s]);
    end
    
    if ~isempty(CLim)
        set(ax, 'CLim', CLim);
    end

    
    
      