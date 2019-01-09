global state
global gh


importRepositories({'Cellbase', 'CSHL', 'mountainsort_matlab_wrapper'});



% % The properties we've been using in the figures
set(0,'defaultLineLineWidth',1);   % set the default line width to lw
set(0,'defaultLineMarkerSize',6); % set the default line marker size to msz



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