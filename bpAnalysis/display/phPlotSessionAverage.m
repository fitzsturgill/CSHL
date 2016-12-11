function varargout = phPlotSessionAverage(analysis, type, outcome, fCh, varargin)

% [ax, hl] = phPlotSessionAverage(analysis, type, outcome, fCh, varargin);

    
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
    hl = zeros(size(type));
    for counter = 1:length(type)
        thisType = type(counter);
        thisOutcome = outcome(counter);
        thisLinespec = linespec{counter};
        avgX = analysis.Photometry.data(thisType, thisOutcome).x(:,1); %x shouldn't really have 2 dimensions, i.e. it is the same for both channels
        avg = analysis.Photometry.data(thisType, thisOutcome).avg(:,fCh);
        avgSEM = analysis.Photometry.data(thisType, thisOutcome).sem(:,fCh);
%         avgX = decimate(avgX, downsample);
%         avg = decimate(avg, downsample);
%         avgSEM = decimate(avgSEM, downsample); % this is wrong
        try
            thisHl = boundedline(avgX, avg, avgSEM, thisLinespec, ax, 'alpha');        
        catch
            thisHl = plot(NaN, NaN); % kludge in case there is only one of that trial type
        end
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
    
        
            