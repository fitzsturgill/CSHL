function saveAllFigures()
% Saves all open figures to several formats in the current directory.
% Accepts a prefix to the filename (figSaveName), to which is appended 
% a) the date and time, and b) 
%
% Modified version of XXX.
% Alex Vaughan, 2015

figHandles   = get (0,'Children');   % locate fall open figure handles
numFigs = size(figHandles,1);
this_now = datestr(now);

if numFigs == 0,
   fprintf('saveAllFigures :: no figures found.') 
end

fprintf('Saving %02g open figures to .pdf, .png, and .fig\n',numFigs)
for f = 1:numFigs
    
    this_gcf = figHandles(f);
    figure(this_gcf)
    figSaveName = sprintf('%s :: %s',get(this_gcf,'Name'),datestr(this_now,'YYYY-mm-DD_HH:MM'));
    figSaveName = strrep(figSaveName,'/','_');
    figSaveName = strrep(figSaveName,'\','_');

    set(this_gcf, 'PaperOrientation', 'portrait');
    if strfind(get(this_gcf,'Name'),'Behavioral Response')
        set(this_gcf, 'PaperSize', [10 15]);
    else
        set(this_gcf, 'PaperSize', [10 8.5]);
    end
    
    fillPage(this_gcf);
    %set(this_gcf,'PaperPositionMode','auto')
    
    fprintf('Saving :: %s\n',figSaveName)
    
    
    
%   imageData = export_fig
%   [imageData, alpha] = export_fig
%   export_fig filename
%   export_fig filename -format1 -format2
%   export_fig ... -nocrop
%   export_fig ... -transparent
%   export_fig ... -native
%   export_fig ... -m<val>
%   export_fig ... -r<val>
%   export_fig ... -a<val>
%   export_fig ... -q<val>
%   export_fig ... -p<val>
%   export_fig ... -d<gs_option>
%   export_fig ... -depsc
%   export_fig ... -<renderer>
%   export_fig ... -<colorspace>
%   export_fig ... -append
%   export_fig ... -bookmark
%   export_fig ... -clipboard
%   export_fig ... -update
%   export_fig ... -nofontswap
%   export_fig(..., handle)

% % 8.24.2015 - can't get ghostscript running, so giving up on this for now.
%     export_fig([figSaveName '.eps'],'-deps','-transparent')
%     export_fig([figSaveName '.png'],'-dpng','-transparent')
%     export_fig([figSaveName '.pdf'],'-dpdf','-transparent')
%     
    print('-deps','-painters', '-r600',[figSaveName '.eps']);
    print('-dpng','-painters', '-r600',[figSaveName '.png']);
    try
        print('-dpdf','-painters', '-r600',[figSaveName '.pdf']);
    catch
        fprintf('Failed to print pdf.\n')
    end
  
    savefig(this_gcf,[figSaveName '.fig'])
    
    
end
fprintf('... done\n')

end

function oldSettings = fillPage(h, varargin)
%FILLPAGE set figure h to fill a page when printed/exported
%   FILLPAGE will adjust printing/exporting related properties so the
%   figure (h) will take up a full page when printed/exported. Optional
%   arguments (varargin) allow you to adjust the margins and page size used.
%
%	OLDSETTINGS = FILLPAGE(fig) returns a structure containing the names of
%	page-related properties that **may** have been affected by FILLPAGE, and
%	the original values of those properties (before any FILLPAGE changes)
%   The figure settings can be restored by 
%      set(fig, oldsettings);
%
%   FILLPAGE also accepts up to 2 optional parameter-value pairs:
%       parameter      value/meaning
%       ----------     --------------------------------------
%        'margins'     1x4 numeric array specifying the 
%                      margins to use on the left, right, top, 
%                      and bottom, respectively. Units used are the 
%                      figure's current PaperUnits.
% 
%                      If not specified, a margin of [.25 .25 .5 .5] 
%                      is used.
%
%        'papersize'   a 1x2 numeric array specifying the 
%                      width and height of the page, respectively. 
%                      Units used are the figure's current PaperUnits.
%                                  OR 
%                      a string specifying the paper type to use 
%                        (see the figure's PaperType property for 
%                        valid values)
%
%                      If not specified, the figure's current PaperSize 
%                      is used.
%
%   Example usage:
%     % fills page based on current settings
%        figure;
%        surf(peaks);
%        oldSettings = fillpage(gcf); 
%     % fills page leaving a margin of .5 around all sides
%        figure;
%        surf(peaks);
%        oldSettings = fillpage(gcf, 'margins', [.5 .5 .5 .5]);
%     % fills page using A4 paper size
%        figure;
%        surf(peaks);
%        oldSettings = fillpage(gcf, 'papersize', 'A4');
%     % fills page using a paper size of 4x6
%        figure;
%        surf(peaks);
%        oldSettings = fillpage(gcf, 'papersize', [4 6]);
%     % fills page using USLetter paper size with margin of 1 all around
%        figure;
%        surf(peaks);
%        oldSettings = fillpage(gcf, 'papersize', 'usletter', 'margins', [1 1 1 1]);
%
% See also PRINT 

%   Copyright 2007 The MathWorks, Inc.

% set up some meaningful aliases for index values
IDX_PAPERTYPE = 1;
IDX_PAPERSIZE = 2;
IDX_PAPERPOSITION = 3;
IDX_PAPERPOSITIONMODE = 4;

if ishandle(h) 
    setprops = {'PaperType', 'PaperSize', 'PaperPosition', 'PaperPositionMode'};
    try 
        settings = get(h, setprops);
    catch 
        error('fillPage unable to get figure properties');
    end
    
    oldSettings = cell2struct(settings, setprops, 2); % to return to caller
    paramSettings = processArgs(varargin{:});
    if isstruct(paramSettings) 
        % process papersize/type
        if isfield(paramSettings, 'papertype')
            if ~isempty(paramSettings.papertype)
                settings{IDX_PAPERTYPE} = paramSettings.papertype;
                % surround with try/catch in case setting is invalid
                try
                   set(h, 'PaperType', paramSettings.papertype);
                catch
                    error('fillPage: invalid PaperType ''%s'' requested ', ...
                        paramSettings.papertype );
                end
                if isfield(paramSettings, 'papersize') && ...
                    ~isempty(paramSettings.papersize) 
                    settings{IDX_PAPERSIZE} = paramSettings.papersize;
                else
                    % since changing papertype might change papersize...
                    % re-read value
                    settings{IDX_PAPERSIZE} = get(h, 'PaperSize'); 
                end
            end
        end
        
        % process margins
        if isfield(paramSettings, 'margins') 
            margins = paramSettings.margins;
            papersize = settings{IDX_PAPERSIZE};
            settings{IDX_PAPERPOSITIONMODE} = 'manual';
            settings{IDX_PAPERPOSITION} = [ ...
                 margins.left ... % X = leftMargin
                 margins.bottom ... % Y = bottomMargin
                 papersize(1) - (margins.right + margins.left) ... % W =  paperwidth - rightMargin - leftMargin
                 papersize(2) - (margins.top + margins.bottom) ... % H = paperheight - topMargin - bottomMargin
                 ];
        end
    end
    set(h, setprops, settings);
else
    error('fillPage requires a handle to a figure');
end

end

% helper function to deal with parsing command line arguments
function param_settings = processArgs(varargin) 
    margins.left   = .25;
    margins.right  = .25;
    margins.top    = .5;
    margins.bottom = .5;
    papertype = []; 
    papersize = []; 

    for i = 1 : 2 : length(varargin)-1 
        param_arg = varargin{i};
        
        % if caller specified margins to use
        if ischar(param_arg) && strcmpi(param_arg, 'margins')
           val_arg = varargin{i+1};
           if isnumeric(val_arg) && length(val_arg) == 4
               margins.left   = val_arg(1);
               margins.right  = val_arg(2);
               margins.top    = val_arg(3);
               margins.bottom = val_arg(4);
           else 
               warning('fillpage:InvalidMargin', ...
                   'fillpage ignoring invalid margin setting; using default');
           end
        end

        % if caller specified papersize to use
        if ischar(param_arg) && strcmpi(param_arg, 'papersize')
           val_arg = varargin{i+1};
           if ischar(val_arg) 
              papertype = val_arg;
            elseif isnumeric(val_arg) && length(val_arg) == 2
              papertype = '<custom>';
              papersize = val_arg;
           else
               warning('fillpage:InvalidPapersize', ...
                   'fillpage ignoring invalid papersize setting; using figure''s current papersize');
           end
        end
    end
    
    param_settings.margins = margins;
    param_settings.papertype = papertype;
    param_settings.papersize = papersize;
end