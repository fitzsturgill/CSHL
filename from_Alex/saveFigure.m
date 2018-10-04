%saveFigure
function saveFigure(figSaveName)

set(gcf, 'PaperOrientation', 'portrait');
set(gcf, 'PaperSize', [22 22]);
set(gcf,'PaperPositionMode','auto')

if nargin <1,
    print('-dpdf','-painters', '-r600',input(sprintf('\nSaving in ... %s.\n What name do you want to save as?\n -> ',pwd),'s'));
else
    print('-dpdf','-painters', '-r600',figSaveName);
end