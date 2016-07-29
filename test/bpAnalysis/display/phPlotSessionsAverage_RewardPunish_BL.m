function varargout = phPlotSessionsAverage_RewardPunish_BL(sessions, typeMatrix, outcomeMatrix, fCh, varargin)

% [ax, hl] = phPlotSessionAverage_RewardPunish_BL(sessions, typeMatrix, outcomeMatrix, fCh, varargin)
% sessions- structure vector containing indv. sessions
% typeMatrix- first row, types to plot, second row, matching types by which
% to baseline each element in the first row
% outcomeMatrix- ditto for outcomes

    
    %% default values
    linespec = {'k', 'r', 'b', 'g'};
%     downsample = 1; % decimate by a factor of 10 (reduces size of figure file)
% NOTE! downsampling of standard error trace is not valid so from now on
% downsample/decimate within processAnalysis_Photometry
    ax = [];
    fig = gcf;
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
%             case 'downsample'
%                 downsample = val;
            case 'ax'
                ax = val;
            case 'fig'
                fig = val;
            otherwise
        end
        counter=counter+2;
    end    
    %%   
    
    if isempty(ax) %make a new axes unless one is provided
        ax=axes(...
        'Parent', fig...
        );
    else
        axes(ax);
    end
    
    hl = zeros(1, size(typeMatrix, 2)); % second row of pair indicates from where to baseline
    for counter = 1:size(typeMatrix, 2)
        thisType = typeMatrix(1, counter);
        thisTypeBL = typeMatrix(2, counter);
        thisOutcome = outcomeMatrix(1, counter);
        thisOutcomeBL = outcomeMatrix(2, counter);
        thisLinespec = linespec{counter};
        % just use the xData from the first session, they should all be the
        % same
        avgX = sessions(1).analysis.Photometry.data(thisType, thisOutcome).x(:,1); %x shouldn't really have 2 dimensions, i.e. it is the same for both channels
        data = [];
        blData = [];
        for j = 1:length(sessions)
            data = [data; sessions(j).analysis.Photometry.data(thisType, thisOutcome).data(:,:,fCh)];
            blData = [blData; sessions(j).analysis.Photometry.data(thisTypeBL, thisOutcomeBL).data(:,:,fCh)];
        end
%         avg = nanmean(data);
%         blAvg = repmat(nanmean(blData));
        data = data - repmat(nanmean(blData), size(data, 1), 1); % baseline
        avg = nanmean(data);
        avgSEM = std(data, 'omitnan') ./ sqrt(sum(~isnan(data), 1));
        thisHl = boundedline(avgX, avg, avgSEM, thisLinespec, ax, 'alpha');        
        hl(counter) = thisHl;
    end
        
    if nargout >=1
        varargout{1} = ax;
    end
    
    if nargout == 2
        varargout{2} = hl;
    end
        
        
        
        
        
        % [ax, 
% nargchk(0, 2, nargout);
% 
% if nargout >= 1
%     varargout{1} = hl;
% end
% 
% if nargout == 2
%     varargout{2} = hp;
% end
    
        
            