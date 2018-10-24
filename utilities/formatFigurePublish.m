function formatFigurePublish(varargin)
    defaults = {...
        'aspect', [2 1];... % aspect ratio w x h (i believe)
        'scaleFactor', 1;...
        'size', []
        };

    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    
    
    fontSize = 8;

    paperUnits = 'inches';
    screenUnits = 'inches';
    if isempty(s.size)
        screenPosition = [0.25 0.25 1 * scaleFactor * aspect(2) 1 * scaleFactor * aspect(1)]; % left bottom width height
    else
        screenPosition = [0.25 0.25 s.size];
    end
    
        
    paperPosition = [0 0 screenPosition(3:4)];
    % in points 8.5 x 11 is 595 x 770

    axs = findobj(gcf, 'Type', 'axes');  
    set(axs, 'TickDir', 'out', 'LineWidth', 1, 'FontName', 'Calibri', 'Box', 'off', 'FontSize', fontSize);
    set(gcf, 'PaperUnits', paperUnits, 'Units', screenUnits);
    set(gcf, 'Position', screenPosition, 'PaperPosition', paperPosition, 'Color', [1 1 1]);
    
    
    
%     https://blogs.mathworks.com/loren/2007/12/11/making-pretty-graphs/

%%
% https://dgleich.wordpress.com/2013/06/04/creating-high-quality-graphics-in-matlab-for-papers-and-presentations/
% % The new defaults will not take effect if there are any open figures. To
% % use them, we close all figures, and then repeat the first example.
% close all;
% 
% % The properties we've been using in the figures
% set(0,'defaultLineLineWidth',1.5);   % set the default line width to lw
% set(0,'defaultLineMarkerSize',8); % set the default line marker size to msz
% 
% % Set the default Size for display
% defpos = get(0,'defaultFigurePosition');
% set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);
% 
% % Set the defaults for saving/printing to a file
% set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
% set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
% defsize = get(gcf, 'PaperSize');
% left = (defsize(1)- width)/2;
% bottom = (defsize(2)- height)/2;
% defsize = [left, bottom, width, height];
% set(0, 'defaultFigurePaperPosition', defsize);
% 
% % Now we repeat the first example but do not need to include anything
% % special beyond manually specifying the tick marks.
% figure(1); clf;
% plot(dmn,f(dmn),'b-',dmn,g(dmn),'r--',xeq,f(xeq),'g*');
% xlim([-pi pi]);
% legend('f(x)', 'g(x)', 'f(x)=g(x)', 'Location', 'SouthEast');
% xlabel('x');
% title('Automatic Example Figure');
% set(gca,'XTick',-3:3); %<- Still need to manually specific tick marks
% set(gca,'YTick',0:10); %<- Still need to manually specific tick marks